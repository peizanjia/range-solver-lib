#include "rstructural_linear_model.h"
#include "rsolverstress.h"
#include "rstructural_control_force_map.h"
#include "rstructural_cable_model_adapter.h"
#include "rstructural_static_truss_model.h"
#include "rmodal_damping.h"

#include <QJsonArray>
#include <algorithm>
#include <cmath>
#include <limits>
#include <rml_model.h>

namespace {
void require(bool ok, const char *message)
{
    if (!ok) throw RError(RError::Type::InvalidInput,R_ERROR_REF,"Stage 9: %s",message);
}
QJsonArray vectorJson(const RRVector &v)
{
    QJsonArray a;
    for (double x : v) { require(std::isfinite(x),"Non-finite vector."); a.append(x); }
    return a;
}
QJsonArray matrixJson(const RRMatrix &m)
{
    QJsonArray a;
    for (uint i=0;i<m.getNRows();++i) a.append(vectorJson(m[i]));
    return a;
}
RRMatrix dense(const RSparseMatrix &s)
{
    RRMatrix m(s.getNRows(),s.getNRows(),0.0);
    for (uint i=0;i<m.getNRows();++i)
        for (uint j=0;j<m.getNRows();++j) m[i][j]=s.findValue(i,j);
    return m;
}
//! Equilibrated SPD factorization; unit-dependent absolute pivots are avoided.
RRVector spdSolve(const RRMatrix &a, const RRVector &b)
{
    const uint n=a.getNRows();
    RRVector scale(n),y(n,0.0),x(n,0.0);
    RRMatrix l(n,n,0.0);
    for (uint i=0;i<n;++i) {
        require(std::isfinite(a[i][i]) && a[i][i]>0.0,"Mass/stiffness has a nonpositive diagonal or mechanism.");
        scale[i]=std::sqrt(a[i][i]);
    }
    for (uint i=0;i<n;++i) for (uint j=0;j<=i;++j) {
        const double denominator=scale[i]*scale[j];
        const double value=a[i][j]/denominator,other=a[j][i]/denominator;
        require(std::isfinite(value) && std::isfinite(other) && std::fabs(value-other)<=1e-10*std::max({1.0,std::fabs(value),std::fabs(other)}),"Mass/stiffness is not finite symmetric.");
        double pivot=value;
        for (uint k=0;k<j;++k) pivot-=l[i][k]*l[j][k];
        if (i==j) { require(std::isfinite(pivot) && pivot>1e-12,"Mass/stiffness is singular, unstable or numerically unresolved."); l[i][j]=std::sqrt(pivot); }
        else l[i][j]=pivot/l[j][j];
    }
    for (uint i=0;i<n;++i) {
        double value=b[i]/scale[i];
        for (uint j=0;j<i;++j) value-=l[i][j]*y[j];
        y[i]=value/l[i][i];
    }
    for (uint row=n;row>0;--row) {
        const uint i=row-1; double value=y[i];
        for (uint j=i+1;j<n;++j) value-=l[j][i]*x[j];
        x[i]=value/l[i][i];
    }
    for (uint i=0;i<n;++i) { x[i]/=scale[i]; require(std::isfinite(x[i]),"Equilibrium solve overflowed."); }
    return x;
}
void constantCondition(const RCondition &condition)
{
    for (uint i=0;i<condition.size();++i) {
        const auto &c=condition.getComponent(i);
        if (!c.getEnabled() || c.empty()) continue;
        const double value=c.getValue(0);
        require(std::isfinite(value),"Boundary/environment conditions must be finite.");
        for (uint j=1;j<c.size();++j)
            require(std::isfinite(c.getValue(j)) && c.getValue(j)==value,"Time-varying loads/supports cannot be frozen silently.");
        if (c.getType()==R_VARIABLE_TEMPERATURE) require(value==293.15,"Thermal working points are outside this exporter.");
    }
}
void validate(const RModel &model)
{
    const auto &p=model.getProblemSetup();
    const auto &interaction=p.getStructuralInteractionSetup();
    const auto mask=model.getProblemTaskTree().getProblemTypeMask();
    require(mask==R_PROBLEM_STRESS || mask==R_PROBLEM_STRESS_MODAL,"Only a structural stress/modal task may be exported.");
    require(!p.getRestart() && !p.getStructuralDynamicsRestartData().available,"Restart states need a separate operating-point binding.");
    require(!p.getStructuralPiezoSetup().enabled,"Dynamic electric/piezo state export is outside this exporter.");
    require(interaction.getLinearConstraints().empty() && interaction.getCouplings().empty() && interaction.getPrescribedMotions().empty(),"MPC/coupling/moving supports need a descriptor constraint formulation.");
    require(!interaction.getActuators().empty() && !interaction.getSensors().empty(),"Export needs explicit input and output catalogs.");
    require(interaction.getActuators().size()<=64 && interaction.getSensors().size()<=64,"Catalog exceeds 64 input/output channels.");
    require(model.getNNodes()>0 && model.getNNodes()<=4096,"Mesh exceeds the bounded dense-export preflight.");
    for (uint e=0;e<model.getNElements();++e) {
        const auto type=model.getElement(e).getType();
        require(type==R_ELEMENT_POINT || type==R_ELEMENT_TRUSS1 || type==R_ELEMENT_BEAM1,"Exporter supports POINT/TRUSS1/BEAM1 only.");
    }
    for (uint g=0;g<model.getNElementGroups();++g) {
        const auto *group=model.getElementGroupPtr(g);
        for (uint b=0;b<group->getNBoundaryConditions();++b) constantCondition(group->getBoundaryCondition(b));
        for (uint e=0;e<group->getNEnvironmentConditions();++e) constantCondition(group->getEnvironmentCondition(e));
        for (uint i=0;i<group->getNInitialConditions();++i) constantCondition(group->getInitialCondition(i));
    }
    const uint temp=model.findVariable(R_VARIABLE_TEMPERATURE);
    if (temp!=RConstants::eod) {
        const auto &v=model.getVariable(temp);
        for (uint i=0;i<v.getNVectors();++i) for (uint j=0;j<v.getNValues();++j)
            require(v.getValue(i,j)==293.15,"Accepted temperature fields must remain at 293.15 K.");
    }
    for (const auto &s:interaction.getSensors()) require(s.aggregation==RStructuralSensorAggregation::Point,"Sensor aggregation must be Point.");
    p.getStructuralDynamicsSetup().validate();
    const auto report=interaction.validate(model);
    require(!report.hasErrors(),"Invalid physical input/output catalog.");
}

class Assembler final : public RSolverStress
{
public:
    Assembler(RModel &model,RSolverSharedData &data):RSolverStress(&model,QString(),QString(),data,false) {}
    QJsonObject extract(const RModel &original,const RStructuralStaticTrussModelData *cableData,
                        const RStructuralCableEquilibriumResult *equilibrium)
    {
        firstRun=true; taskIteration=0; scales=RScales();
        preparingInitialForce=true; // Assemble once; no t0/t1 transient force prepass.
        findComputableElements(problemType); findIncludableElements(); findInwardElements(); updateLocalRotations();
        recover(); generateLocalConstraints(); generateEquationBook();
        require(equationBook.getNEnabled()>0 && equationBook.getNEnabled()<=RStructuralLinearModel::MaxEquations,"Free equation count must be 1..256.");
        prepare();
        const uint n=A.getNRows();
        RRVector q(n,0.0);
        double residual=0.0;
        if (cableData) {
            A.clear(); A.setNRows(n);
            M.clear(); M.setNRows(n);
            const auto &setup=original.getProblemSetup().getStructuralCableSetup();
            RRVector nodeMass(original.getNNodes(),0.0);
            for (const auto &member:cableData->members) {
                const RStructuralCableMember *config=nullptr;
                for (const auto &candidate:setup.members) if (candidate.elementId==member.elementId) config=&candidate;
                require(config!=nullptr,"Cable member binding disappeared.");
                auto p0=cableData->nodes[member.node1],p1=cableData->nodes[member.node2];
                for (uint axis=0;axis<3;++axis) {
                    p0[axis]+=equilibrium->displacements[3*member.node1+axis];
                    p1[axis]+=equilibrium->displacements[3*member.node2+axis];
                }
                const auto kernel=RStructuralCableKernel::evaluate(p0,p1,config->referenceLength,config->initialTension,member.area*member.youngModulus);
                require(kernel.taut && kernel.trialTension>1e-10*std::max({1.0,member.area*member.youngModulus,std::fabs(config->initialTension)}),"Slack/taut cutoff or unresolved branch cannot be smoothly linearized.");
                const uint nodes[2]={member.node1,member.node2};
                for (uint ni=0;ni<2;++ni) for (uint nj=0;nj<2;++nj)
                    for (uint di=0;di<3;++di) for (uint dj=0;dj<3;++dj) {
                        uint row=0,col=0;
                        if (!equationBook.getValue(nodes[ni],RDofType(di),row) || !equationBook.getValue(nodes[nj],RDofType(dj),col)) continue;
                        const auto ai=axis(nodes[ni],RDofType(di)),aj=axis(nodes[nj],RDofType(dj));
                        double value=0.0;
                        for (uint i=0;i<3;++i) for (uint j=0;j<3;++j) value+=ai[i]*kernel.tangent[3*ni+i][3*nj+j]*aj[j];
                        A.addValue(row,col,value);
                    }
                double density=0.0;
                for (uint lineId=0;lineId<original.getNLines();++lineId) {
                    const auto &line=original.getLine(lineId);
                    for (uint e=0;e<line.size();++e) if (line.get(e)==member.elementId) {
                        const auto &material=line.getMaterial();
                        const uint position=material.findPosition(RMaterialProperty::Density);
                        require(position<material.size(),"Cable export requires positive mass density.");
                        density=material.get(position).get(293.15);
                    }
                }
                require(std::isfinite(density) && density>0.0,"Cable export requires positive mass density.");
                const double halfMass=0.5*density*member.area*member.referenceLength;
                require(std::isfinite(halfMass) && halfMass>0.0,"Cable reference-geometry mass is invalid.");
                nodeMass[member.node1]+=halfMass; nodeMass[member.node2]+=halfMass;
            }
            for (const auto &dof:orderedPhysicalDofs()) {
                uint row=0; equationBook.getValue(dof.nodeId,dof.dof,row);
                const auto direction=axis(dof.nodeId,dof.dof);
                for (uint component=0;component<3;++component) q[row]+=direction[component]*equilibrium->displacements[3*dof.nodeId+component];
                M.addValue(row,row,nodeMass[dof.nodeId]);
            }
            residual=equilibrium->residualNorm;
        }
        const auto kd=dense(A),md=dense(M);
        spdSolve(md,RRVector(n,0.0));
        const auto solved=spdSolve(kd,b);
        if (!cableData) {
            q=solved;
            RRVector r(n,0.0); double scale=0.0;
            for (uint i=0;i<n;++i) { r[i]=-b[i]; scale=std::max(scale,std::fabs(b[i])); for (uint j=0;j<n;++j) { r[i]+=kd[i][j]*q[j]; scale=std::max(scale,std::fabs(kd[i][j]*q[j])); } }
            residual=RRVector::euclideanNorm(r);
            require(residual<=1e-9*std::max(scale,std::numeric_limits<double>::min())*std::sqrt(double(n)),"Fresh equilibrium residual certificate failed.");
        }
        prepareControlMapping();
        const auto &catalog=original.getProblemSetup().getStructuralInteractionSetup();
        const uint m=uint(catalog.getActuators().size()),p=uint(controlSensors.size());
        RRMatrix input(n,m,0.0),sq(p,n,0.0),sv(p,n,0.0),sa(p,n,0.0),direct(p,m,0.0);
        for (uint j=0;j<m;++j) for (const auto &weight:controlForceMap->getChannels()[j].reducedWeights) input[weight.first][j]=weight.second;
        RRVector output(p,0.0);
        QJsonArray outputs,inputs,dofs;
        for (uint i=0;i<p;++i) {
            const auto &s=controlSensors[i];
            auto &row=s.quantity==RStructuralSensorQuantity::Displacement ? sq[i] : s.quantity==RStructuralSensorQuantity::Velocity ? sv[i] : sa[i];
            for (const auto &w:s.weights) row[w.first]=w.second;
            if (s.quantity==RStructuralSensorQuantity::Displacement) { output[i]=s.displacementOffset; for (const auto &w:s.weights) output[i]+=w.second*q[w.first]; }
            const auto &target=catalog.getSensors()[i].target;
            const QString quantity=s.quantity==RStructuralSensorQuantity::Displacement ? "Displacement" : s.quantity==RStructuralSensorQuantity::Velocity ? "Velocity" : "Acceleration";
            outputs.append(QJsonObject{{"id",s.id},{"unit",s.unit},{"quantity",quantity},{"nodePosition",int(target.nodeId)},{"nodePhysicalId",QString::number(original.getNodePhysicalId(target.nodeId))}});
        }
        for (const auto &a:catalog.getActuators()) inputs.append(QJsonObject{{"id",a.id},{"unit",a.unit},{"nodePosition",int(a.target.nodeId)},{"nodePhysicalId",QString::number(original.getNodePhysicalId(a.target.nodeId))}});
        const char *names[]={"UX","UY","UZ","RX","RY","RZ"};
        for (const auto &d:orderedPhysicalDofs()) {
            const auto direction=axis(d.nodeId,d.dof);
            const auto &rotation=uint(d.dof)<3 ? localRotations[d.nodeId] : rotationLocalRotations[d.nodeId];
            dofs.append(QJsonObject{{"nodePosition",int(d.nodeId)},{"nodePhysicalId",QString::number(original.getNodePhysicalId(d.nodeId))},{"dof",names[uint(d.dof)]},{"unit",uint(d.dof)<3 ? "m" : "rad"},{"frame",rotation.isActive() ? "constraint-local" : "global"},{"globalDirection",QJsonArray{direction[0],direction[1],direction[2]}}});
        }
        const auto &dyn=original.getProblemSetup().getStructuralDynamicsSetup();
        RRMatrix damping(n,n,0.0);
        for (uint i=0;i<n;++i) for (uint j=0;j<n;++j) damping[i][j]=dyn.getRayleighAlpha()*md[i][j]+dyn.getRayleighBeta()*kd[i][j];
        if (!dyn.getModalDampingRatios().empty()) {
            const auto modal=RModalDamping::synthesize(M,A,dyn.getModalDampingRatios(),original.getMatrixSolverConf(RMatrixSolverConf::CG));
            for (uint i=0;i<n;++i) for (uint j=0;j<n;++j) damping[i][j]+=modal.damping.findValue(i,j);
        }
        return {{"schema","range-fea-linear-model"},{"version",1},{"M",matrixJson(md)},{"C",matrixJson(damping)},{"K",matrixJson(kd)},{"B",matrixJson(input)},
                {"Sq",matrixJson(sq)},{"Sv",matrixJson(sv)},{"Sa",matrixJson(sa)},{"D",matrixJson(direct)},
                {"operating_point",QJsonObject{{"q",vectorJson(q)},{"v",vectorJson(RRVector(n,0.0))},{"a",vectorJson(RRVector(n,0.0))},{"force",vectorJson(b)},{"outputs",vectorJson(output)},{"residual",residual}}},
                {"inputs",inputs},{"outputs",outputs},{"dofs",dofs},{"source",QJsonObject{{"kind",cableData ? "taut-cable-tangent-FE" : "linear-elastic-FE"},{"units","SI"},{"coordinateMeaning","perturbations about static equilibrium"},{"massConvention",cableData ? "reference-geometry lumped" : original.getProblemSetup().getMassMatrixType()==RMassMatrixType::Lumped ? "lumped" : "consistent"}}}};
    }
private:
    RR3Vector axis(uint node,RDofType dof) const {
        const uint component=uint(dof)%3;
        const auto &rotation=uint(dof)<3 ? localRotations[node] : rotationLocalRotations[node];
        RR3Vector direction(0.0,0.0,0.0);
        if (!rotation.isActive()) direction[component]=1.0;
        else for (uint i=0;i<3;++i) direction[i]=rotation.getR()[i][component];
        return direction;
    }
};
}

QJsonObject RStructuralLinearModel::build(const RModel &model)
{
    validate(model);
    RModel copy(model);
    const bool cable=model.getProblemSetup().getStructuralCableSetup().enabled;
    RStructuralStaticTrussModelData data;
    RStructuralCableEquilibriumResult equilibrium;
    if (cable) {
        require(model.getProblemSetup().getStructuralDynamicsSetup().getRayleighAlpha()==0.0 && model.getProblemSetup().getStructuralDynamicsSetup().getRayleighBeta()==0.0 && model.getProblemSetup().getStructuralDynamicsSetup().getModalDampingRatios().empty(),"Cable dynamics has no Rayleigh/modal damping path.");
        RModel staticModel(model);
        staticModel.getProblemSetup().getStructuralDynamicsSetup().setEnabled(false);
        staticModel.getTimeSolver().setEnabled(false);
        staticModel.getProblemSetup().getStructuralControlSetup()=RStructuralControlSetup();
        staticModel.getProblemSetup().getStructuralInteractionSetup().getActuators().clear();
        staticModel.getProblemSetup().getStructuralInteractionSetup().getSensors().clear();
        data=RStructuralStaticTrussModel::build(staticModel);
        for (uint i=0;i<data.prescribedMask.size();++i) if (data.prescribedMask[i]) require(data.prescribedDisplacements[i]==0.0,"Cable export requires fixed zero-valued supports.");
        equilibrium=RStructuralCableModelAdapter::solve(staticModel);
        copy.getProblemSetup().getStructuralCableSetup().enabled=false;
        copy.getProblemSetup().setMassMatrixType(RMassMatrixType::Lumped);
    }
    // The transient stress assembly supplies physical M/K and the correct
    // static held-column RHS and moment loads. Do not invoke its solve or
    // step lifecycle; legacy time integration and modal prestress geometry
    // are deliberately excluded from this reference-geometry linearization.
    copy.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    copy.getProblemSetup().getStructuralDynamicsSetup().setEnabled(true);
    copy.getTimeSolver().setEnabled(false);
    RSolverSharedData shared;
    Assembler assembler(copy,shared);
    return assembler.extract(model,cable ? &data : nullptr,cable ? &equilibrium : nullptr);
}
