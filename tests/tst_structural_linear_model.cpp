#include <QtTest>

#include <QJsonArray>
#include <QJsonObject>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <rbl_error.h>
#include <rml_model.h>
#include <rml_structural_interaction.h>
#include <rml_variable.h>

#include "rstructural_linear_model.h"
#include "stage9_test_models.h"

namespace
{
using Stage9TestModels::Area;
using Stage9TestModels::CableTension;
using Stage9TestModels::Density;
using Stage9TestModels::Length;
using Stage9TestModels::YoungModulus;

bool near(double actual,double expected,double relativeTolerance = 1.0e-9,
          double absoluteTolerance = 1.0e-11)
{
    return std::isfinite(actual) && std::isfinite(expected) &&
        std::fabs(actual-expected) <= absoluteTolerance+
            relativeTolerance*std::max(std::fabs(actual),std::fabs(expected));
}

double number(const QJsonValue &value)
{
    return value.toDouble(std::numeric_limits<double>::quiet_NaN());
}

double matrixAt(const QJsonObject &object,const QString &name,int row,int column)
{
    const QJsonArray matrix = object.value(name).toArray();
    if (row < 0 || row >= matrix.size()) return std::numeric_limits<double>::quiet_NaN();
    const QJsonArray values = matrix.at(row).toArray();
    if (column < 0 || column >= values.size()) return std::numeric_limits<double>::quiet_NaN();
    return number(values.at(column));
}

double vectorAt(const QJsonObject &object,const QString &name,int position)
{
    const QJsonArray values = object.value(name).toArray();
    if (position < 0 || position >= values.size()) return std::numeric_limits<double>::quiet_NaN();
    return number(values.at(position));
}

QJsonObject operatingPoint(const QJsonObject &model)
{
    return model.value(QStringLiteral("operating_point")).toObject();
}

int outputPosition(const QJsonObject &model,const QString &id)
{
    const QJsonArray outputs = model.value(QStringLiteral("outputs")).toArray();
    for (int i=0;i<outputs.size();++i)
    {
        if (outputs.at(i).toObject().value(QStringLiteral("id")).toString() == id) return i;
    }
    return -1;
}

int dofPosition(const QJsonObject &model,uint nodePosition,
                const QString &unit,const std::array<double,3> &globalDirection)
{
    const QJsonArray dofs = model.value(QStringLiteral("dofs")).toArray();
    for (int i=0;i<dofs.size();++i)
    {
        const QJsonObject dof = dofs.at(i).toObject();
        if (dof.value(QStringLiteral("nodePosition")).toInt(-1) != int(nodePosition) ||
            dof.value(QStringLiteral("unit")).toString() != unit)
        {
            continue;
        }
        const QJsonArray direction = dof.value(QStringLiteral("globalDirection")).toArray();
        if (direction.size() != 3) continue;
        bool matches = true;
        for (int axis=0;axis<3;++axis)
        {
            matches = matches && near(number(direction.at(axis)),globalDirection[axis],1.0e-8,1.0e-10);
        }
        if (matches) return i;
    }
    return -1;
}

bool sameMatrix(const QJsonObject &a,const QJsonObject &b,const QString &name)
{
    const QJsonArray am = a.value(name).toArray();
    const QJsonArray bm = b.value(name).toArray();
    if (am.size() != bm.size()) return false;
    for (int row=0;row<am.size();++row)
    {
        const QJsonArray ar = am.at(row).toArray();
        const QJsonArray br = bm.at(row).toArray();
        if (ar.size() != br.size()) return false;
        for (int column=0;column<ar.size();++column)
        {
            if (!near(number(ar.at(column)),number(br.at(column)),1.0e-12,1.0e-13)) return false;
        }
    }
    return true;
}

bool isScaledSpd(const QJsonObject &model,const QString &name,int size)
{
    const QJsonArray matrix = model.value(name).toArray();
    if (matrix.size() != size) return false;
    std::vector<double> diagonal(size,0.0);
    std::vector<std::vector<double>> factor(size,std::vector<double>(size,0.0));
    for (int i=0;i<size;++i)
    {
        const QJsonArray row = matrix.at(i).toArray();
        if (row.size() != size) return false;
        diagonal[i] = number(row.at(i));
        if (!std::isfinite(diagonal[i]) || diagonal[i] <= 0.0) return false;
    }
    for (int i=0;i<size;++i)
    {
        for (int j=0;j<=i;++j)
        {
            const double aij = matrixAt(model,name,i,j)/std::sqrt(diagonal[i]*diagonal[j]);
            const double aji = matrixAt(model,name,j,i)/std::sqrt(diagonal[i]*diagonal[j]);
            if (!near(aij,aji,1.0e-9,1.0e-11)) return false;
            double pivot = aij;
            for (int k=0;k<j;++k) pivot -= factor[i][k]*factor[j][k];
            if (i == j)
            {
                if (!std::isfinite(pivot) || pivot <= 1.0e-12) return false;
                factor[i][j] = std::sqrt(pivot);
            }
            else
            {
                factor[i][j] = pivot/factor[j][j];
            }
        }
    }
    return true;
}

void addDisplacementSensor(RModel &model,const QString &id,
                           const RR3Vector &direction)
{
    RStructuralSensor sensor;
    sensor.id = id;
    sensor.target.nodeId = model.getNNodes()-1;
    sensor.quantity = RStructuralSensorQuantity::Displacement;
    sensor.direction = direction;
    sensor.unit = RStructuralInteractionSetup::getSensorUnit(sensor.quantity);
    model.getProblemSetup().getStructuralInteractionSetup().getSensors().push_back(sensor);
}

std::array<double,3> cableResidual(const std::array<double,3> &tipDisplacement)
{
    const double x = Length+tipDisplacement[0];
    const double y = tipDisplacement[1];
    const double z = tipDisplacement[2];
    const double currentLength = std::hypot(std::hypot(x,y),z);
    const double tension = std::max(0.0,CableTension+
        YoungModulus*Area/Length*(currentLength-Length));
    return {{tension*x/currentLength-CableTension,
             tension*y/currentLength,
             tension*z/currentLength}};
}

void addCableAxisChannels(RModel &model)
{
    auto &interaction = model.getProblemSetup().getStructuralInteractionSetup();
    const RVariableType forceTypes[2] = {R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    const RR3Vector directions[2] = {RR3Vector(0.0,1.0,0.0),RR3Vector(0.0,0.0,1.0)};
    const QString inputIds[2] = {QStringLiteral("tip-force-y"),QStringLiteral("tip-force-z")};
    for (int axis=0;axis<2;++axis)
    {
        RStructuralActuator actuator;
        actuator.id = inputIds[axis];
        actuator.target.nodeId = 1;
        actuator.variable = forceTypes[axis];
        actuator.direction = directions[axis];
        actuator.unit = RVariable::getUnits(actuator.variable);
        actuator.minimum = -1.0e6;
        actuator.maximum = 1.0e6;
        interaction.getActuators().push_back(actuator);
    }

    interaction.getSensors().clear();
    const RStructuralSensorQuantity quantities[3] = {
        RStructuralSensorQuantity::Displacement,
        RStructuralSensorQuantity::Velocity,
        RStructuralSensorQuantity::Acceleration
    };
    const char *quantityNames[3] = {"displacement","velocity","acceleration"};
    const RR3Vector axes[3] = {
        RR3Vector(1.0,0.0,0.0),RR3Vector(0.0,1.0,0.0),RR3Vector(0.0,0.0,1.0)
    };
    for (int quantity=0;quantity<3;++quantity)
    {
        for (int axis=0;axis<3;++axis)
        {
            RStructuralSensor sensor;
            sensor.id = QStringLiteral("tip-%1-%2")
                            .arg(QString::fromLatin1(quantityNames[quantity])).arg(axis);
            sensor.target.nodeId = 1;
            sensor.quantity = quantities[quantity];
            sensor.direction = axes[axis];
            sensor.unit = RStructuralInteractionSetup::getSensorUnit(sensor.quantity);
            interaction.getSensors().push_back(sensor);
        }
    }
}
}

class TestStructuralLinearModel : public QObject
{
    Q_OBJECT

private slots:
    void trussExportsNativeSiOperatorsAndPreservesSource();
    void constantTipLoadHasSolvedEquilibriumAndFreshResidual();
    void nonzeroPrescribedSupportAndOldHistoryRemainPhysicalAndDeterministic();
    void frameCantileverMatchesForceAndMomentComplianceWithSpdMass();
    void tautCableExportsPrestressTangentMassAndIndependentJacobian();
    void localFrameDofLabelsAndControlDirectionsResolveGlobally();
    void unsupportedAndUnstableInputsFailClosed();
};

void TestStructuralLinearModel::trussExportsNativeSiOperatorsAndPreservesSource()
{
    RModel model = Stage9TestModels::makeTruss();
    RStructuralDynamicsSetup &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setRayleighAlpha(0.25);
    dynamics.setRayleighBeta(2.0e-3);

    const uint sourceVariables = model.getNVariables();
    const double sourceTipX = model.getNode(1).getX();
    const double sourceTime = model.getTimeSolver().getComputedTime();
    const QJsonObject exported = RStructuralLinearModel::build(model);
    const QJsonObject op = operatingPoint(exported);

    QCOMPARE(exported.value(QStringLiteral("schema")).toString(),QStringLiteral("range-fea-linear-model"));
    QCOMPARE(exported.value(QStringLiteral("version")).toInt(),1);
    QCOMPARE(exported.value(QStringLiteral("source")).toObject()
                 .value(QStringLiteral("units")).toString(),QStringLiteral("SI"));
    QCOMPARE(exported.value(QStringLiteral("source")).toObject()
                 .value(QStringLiteral("massConvention")).toString(),QStringLiteral("consistent"));
    QCOMPARE(int(exported.value(QStringLiteral("dofs")).toArray().size()),1);
    QCOMPARE(int(exported.value(QStringLiteral("inputs")).toArray().size()),1);
    QCOMPARE(int(exported.value(QStringLiteral("outputs")).toArray().size()),3);

    const double expectedMass = Density*Area*Length/3.0;
    const double expectedStiffness = YoungModulus*Area/Length;
    QVERIFY(near(matrixAt(exported,QStringLiteral("M"),0,0),expectedMass));
    QVERIFY(near(matrixAt(exported,QStringLiteral("K"),0,0),expectedStiffness));
    QVERIFY(near(matrixAt(exported,QStringLiteral("C"),0,0),
                   dynamics.getRayleighAlpha()*expectedMass+
                       dynamics.getRayleighBeta()*expectedStiffness));
    QVERIFY(near(matrixAt(exported,QStringLiteral("B"),0,0),1.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sq"),0,0),1.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sq"),1,0),0.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sq"),2,0),0.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sv"),0,0),0.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sv"),1,0),1.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sv"),2,0),0.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sa"),0,0),0.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sa"),1,0),0.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sa"),2,0),1.0));
    for (int row=0;row<3;++row) QVERIFY(near(matrixAt(exported,QStringLiteral("D"),row,0),0.0));

    const QJsonObject dof = exported.value(QStringLiteral("dofs")).toArray().at(0).toObject();
    QCOMPARE(dof.value(QStringLiteral("nodePosition")).toInt(),1);
    QCOMPARE(dof.value(QStringLiteral("nodePhysicalId")).toString(),QStringLiteral("7118"));
    QCOMPARE(dof.value(QStringLiteral("unit")).toString(),QStringLiteral("m"));
    QCOMPARE(dof.value(QStringLiteral("frame")).toString(),QStringLiteral("constraint-local"));
    const QJsonArray direction = dof.value(QStringLiteral("globalDirection")).toArray();
    QVERIFY(near(number(direction.at(0)),1.0));
    QVERIFY(near(number(direction.at(1)),0.0));
    QVERIFY(near(number(direction.at(2)),0.0));

    QCOMPARE(int(op.value(QStringLiteral("q")).toArray().size()),1);
    QCOMPARE(int(op.value(QStringLiteral("v")).toArray().size()),1);
    QCOMPARE(int(op.value(QStringLiteral("a")).toArray().size()),1);
    QVERIFY(near(vectorAt(op,QStringLiteral("q"),0),0.0));
    QVERIFY(near(vectorAt(op,QStringLiteral("v"),0),0.0));
    QVERIFY(near(vectorAt(op,QStringLiteral("a"),0),0.0));
    QVERIFY(near(vectorAt(op,QStringLiteral("force"),0),0.0));
    QVERIFY(number(op.value(QStringLiteral("residual"))) <= 1.0e-12);

    QCOMPARE(model.getNVariables(),sourceVariables);
    QCOMPARE(model.findVariable(R_VARIABLE_DISPLACEMENT),RConstants::eod);
    QVERIFY(near(model.getNode(1).getX(),sourceTipX));
    QVERIFY(near(model.getTimeSolver().getComputedTime(),sourceTime));
    QVERIFY(!model.getTimeSolver().getEnabled());
    QCOMPARE(model.getProblemSetup().getStructuralInteractionSetup().getActuators().size(),std::size_t(1));
    QCOMPARE(model.getProblemSetup().getStructuralInteractionSetup().getSensors().size(),std::size_t(3));
}

void TestStructuralLinearModel::constantTipLoadHasSolvedEquilibriumAndFreshResidual()
{
    RModel model = Stage9TestModels::makeTruss(1,100.0);
    const QJsonObject exported = RStructuralLinearModel::build(model);
    const QJsonObject op = operatingPoint(exported);
    const double q = vectorAt(op,QStringLiteral("q"),0);
    const double force = vectorAt(op,QStringLiteral("force"),0);
    const double freshResidual = matrixAt(exported,QStringLiteral("K"),0,0)*q-force;
    QVERIFY(near(q,0.1,1.0e-8,1.0e-11));
    QVERIFY(near(force,100.0));
    QVERIFY(std::fabs(freshResidual) <= 1.0e-8);
    QVERIFY(std::fabs(number(op.value(QStringLiteral("residual")))) <= 1.0e-8);
    QVERIFY(near(vectorAt(op,QStringLiteral("outputs"),0),0.1,1.0e-8,1.0e-11));
}

void TestStructuralLinearModel::nonzeroPrescribedSupportAndOldHistoryRemainPhysicalAndDeterministic()
{
    RModel base = Stage9TestModels::makeTruss(1,100.0,0.025);
    const QJsonObject baseExport = RStructuralLinearModel::build(base);
    const QJsonObject baseOp = operatingPoint(baseExport);
    QVERIFY(near(vectorAt(baseOp,QStringLiteral("q"),0),0.125,1.0e-8,1.0e-11));
    QVERIFY(near(vectorAt(baseOp,QStringLiteral("outputs"),0),0.125,1.0e-8,1.0e-11));
    QVERIFY(std::fabs(number(baseOp.value(QStringLiteral("residual")))) <= 1.0e-8);

    RModel history = base;
    RTimeSolver &time = history.getTimeSolver();
    time.setEnabled(true);
    time.setTimes({0.0,0.5,1.0});
    time.setCurrentTimeStep(2);
    time.setComputedTime(1.0);
    RVariable oldDisplacement(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE);
    oldDisplacement.resize(3,history.getNNodes());
    oldDisplacement.setValue(0,0,0.4);
    oldDisplacement.setValue(1,0,-0.2);
    oldDisplacement.setValue(2,0,0.1);
    oldDisplacement.setValue(0,1,-0.7);
    oldDisplacement.setValue(1,1,0.3);
    oldDisplacement.setValue(2,1,-0.5);
    history.addVariable(oldDisplacement);
    const QJsonObject historyExport = RStructuralLinearModel::build(history);
    for (const QString &name : {QStringLiteral("M"),QStringLiteral("C"),QStringLiteral("K"),
                                QStringLiteral("B"),QStringLiteral("Sq"),QStringLiteral("Sv"),
                                QStringLiteral("Sa"),QStringLiteral("D")})
    {
        QVERIFY2(sameMatrix(baseExport,historyExport,name),qPrintable(name));
    }
    const QJsonObject historyOp = operatingPoint(historyExport);
    QVERIFY(near(vectorAt(historyOp,QStringLiteral("q"),0),0.125,1.0e-8,1.0e-11));
    QVERIFY(near(vectorAt(historyOp,QStringLiteral("outputs"),0),0.125,1.0e-8,1.0e-11));
    QCOMPARE(history.getNVariables(),1u);
    QCOMPARE(history.findVariable(R_VARIABLE_DISPLACEMENT),0u);
    QVERIFY(near(history.getVariable(0).getValue(0,1),-0.7));
    QVERIFY(near(history.getVariable(0).getValue(1,1),0.3));
    QVERIFY(near(history.getVariable(0).getValue(2,1),-0.5));
    QVERIFY(history.getTimeSolver().getEnabled());
    QCOMPARE(history.getTimeSolver().getCurrentTimeStep(),2u);
    QVERIFY(near(history.getTimeSolver().getComputedTime(),1.0));
}

void TestStructuralLinearModel::frameCantileverMatchesForceAndMomentComplianceWithSpdMass()
{
    constexpr double iz = 5.0e-5;
    constexpr double ei = YoungModulus*iz;
    constexpr double forceY = 0.1;
    RModel forceModel = Stage9TestModels::makeFrameCantilever(forceY,0.0);
    addDisplacementSensor(forceModel,QStringLiteral("tip-displacement-y"),RR3Vector(0.0,1.0,0.0));
    const QJsonObject forceExport = RStructuralLinearModel::build(forceModel);
    const int uy = dofPosition(forceExport,1,QStringLiteral("m"),{{0.0,1.0,0.0}});
    QVERIFY(uy >= 0);
    QVERIFY(near(vectorAt(operatingPoint(forceExport),QStringLiteral("q"),uy),
                 forceY*std::pow(Length,3)/(3.0*ei),1.0e-7,1.0e-10));
    const int ySensor = outputPosition(forceExport,QStringLiteral("tip-displacement-y"));
    QVERIFY(ySensor >= 0);
    QVERIFY(near(vectorAt(operatingPoint(forceExport),QStringLiteral("outputs"),ySensor),
                 forceY*std::pow(Length,3)/(3.0*ei),1.0e-7,1.0e-10));
    QCOMPARE(int(forceExport.value(QStringLiteral("M")).toArray().size()),6);
    QCOMPARE(int(forceExport.value(QStringLiteral("K")).toArray().size()),6);
    QVERIFY(isScaledSpd(forceExport,QStringLiteral("M"),6));

    constexpr double tipMoment = 0.2;
    RModel momentModel = Stage9TestModels::makeFrameCantilever(0.0,tipMoment);
    addDisplacementSensor(momentModel,QStringLiteral("tip-displacement-y"),RR3Vector(0.0,1.0,0.0));
    const QJsonObject momentExport = RStructuralLinearModel::build(momentModel);
    const int momentUy = dofPosition(momentExport,1,QStringLiteral("m"),{{0.0,1.0,0.0}});
    const int rz = dofPosition(momentExport,1,QStringLiteral("rad"),{{0.0,0.0,1.0}});
    QVERIFY(momentUy >= 0);
    QVERIFY(rz >= 0);
    QVERIFY(near(vectorAt(operatingPoint(momentExport),QStringLiteral("q"),momentUy),
                 tipMoment*Length*Length/(2.0*ei),1.0e-7,1.0e-10));
    QVERIFY(near(vectorAt(operatingPoint(momentExport),QStringLiteral("q"),rz),
                 tipMoment*Length/ei,1.0e-7,1.0e-10));
}

void TestStructuralLinearModel::tautCableExportsPrestressTangentMassAndIndependentJacobian()
{
    RModel model = Stage9TestModels::makeTautCable();
    addCableAxisChannels(model);
    const uint sourceVariables = model.getNVariables();
    const QJsonObject exported = RStructuralLinearModel::build(model);
    const QJsonObject op = operatingPoint(exported);
    const int n = int(exported.value(QStringLiteral("dofs")).toArray().size());
    QCOMPARE(n,3);
    QCOMPARE(exported.value(QStringLiteral("source")).toObject()
                 .value(QStringLiteral("massConvention")).toString(),
             QStringLiteral("reference-geometry lumped"));
    QCOMPARE(exported.value(QStringLiteral("source")).toObject()
                 .value(QStringLiteral("kind")).toString(),QStringLiteral("taut-cable-tangent-FE"));
    for (int axis=0;axis<3;++axis)
    {
        QVERIFY(near(matrixAt(exported,QStringLiteral("M"),axis,axis),
                     Density*Area*Length/2.0));
        for (int other=0;other<3;++other)
        {
            const double expected = axis == other
                ? (axis == 0 ? YoungModulus*Area/Length : CableTension/Length)
                : 0.0;
            QVERIFY(near(matrixAt(exported,QStringLiteral("K"),axis,other),expected,1.0e-8,1.0e-10));
            QVERIFY(near(matrixAt(exported,QStringLiteral("B"),axis,other),axis == other ? 1.0 : 0.0));
        }
    }
    for (int row=0;row<9;++row)
    {
        for (int column=0;column<3;++column)
        {
            const bool collocated = row%3 == column;
            QVERIFY(near(matrixAt(exported,QStringLiteral("Sq"),row,column),row<3 && collocated ? 1.0 : 0.0));
            QVERIFY(near(matrixAt(exported,QStringLiteral("Sv"),row,column),row>=3 && row<6 && collocated ? 1.0 : 0.0));
            QVERIFY(near(matrixAt(exported,QStringLiteral("Sa"),row,column),row>=6 && collocated ? 1.0 : 0.0));
            QVERIFY(near(matrixAt(exported,QStringLiteral("D"),row,column),0.0));
        }
    }
    QVERIFY(near(vectorAt(op,QStringLiteral("q"),0),0.0));
    QVERIFY(near(vectorAt(op,QStringLiteral("q"),1),0.0));
    QVERIFY(near(vectorAt(op,QStringLiteral("q"),2),0.0));
    QVERIFY(near(vectorAt(op,QStringLiteral("force"),0),CableTension));
    QVERIFY(near(number(op.value(QStringLiteral("residual"))),0.0,1.0e-8,1.0e-12));

    // Independent central differences of the nonlinear cable residual check
    // each free tip direction at the balanced, undeformed operating point.
    const double step = 1.0e-5;
    for (int column=0;column<3;++column)
    {
        std::array<double,3> plus{{0.0,0.0,0.0}};
        std::array<double,3> minus{{0.0,0.0,0.0}};
        plus[column] = step;
        minus[column] = -step;
        const auto rp = cableResidual(plus);
        const auto rm = cableResidual(minus);
        for (int row=0;row<3;++row)
        {
            const double finiteDifference = (rp[row]-rm[row])/(2.0*step);
            QVERIFY(near(finiteDifference,matrixAt(exported,QStringLiteral("K"),row,column),
                         2.0e-6,1.0e-8));
        }
    }

    QCOMPARE(model.getNVariables(),sourceVariables);
    QVERIFY(model.getProblemSetup().getStructuralCableSetup().enabled);
    QCOMPARE(model.getProblemSetup().getStructuralCableSetup().members.size(),std::size_t(1));
    QVERIFY(near(model.getProblemSetup().getStructuralCableSetup().members[0].initialTension,CableTension));
    QVERIFY(near(model.getNode(1).getX(),Length));
    QCOMPARE(model.getProblemSetup().getStructuralInteractionSetup().getActuators().size(),std::size_t(3));
    QCOMPARE(model.getProblemSetup().getStructuralInteractionSetup().getSensors().size(),std::size_t(9));
}

void TestStructuralLinearModel::localFrameDofLabelsAndControlDirectionsResolveGlobally()
{
    RModel model = Stage9TestModels::makeTruss();
    Stage9TestModels::rotateTrussToY(model);
    const QJsonObject exported = RStructuralLinearModel::build(model);
    QCOMPARE(int(exported.value(QStringLiteral("dofs")).toArray().size()),1);
    const QJsonObject dof = exported.value(QStringLiteral("dofs")).toArray().at(0).toObject();
    QCOMPARE(dof.value(QStringLiteral("nodePhysicalId")).toString(),QStringLiteral("7118"));
    QCOMPARE(dof.value(QStringLiteral("frame")).toString(),QStringLiteral("constraint-local"));
    QVERIFY(near(number(dof.value(QStringLiteral("globalDirection")).toArray().at(0)),0.0));
    QVERIFY(near(number(dof.value(QStringLiteral("globalDirection")).toArray().at(1)),1.0));
    QVERIFY(near(number(dof.value(QStringLiteral("globalDirection")).toArray().at(2)),0.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("K"),0,0),YoungModulus*Area/Length));
    QVERIFY(near(matrixAt(exported,QStringLiteral("B"),0,0),1.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sq"),0,0),1.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sv"),1,0),1.0));
    QVERIFY(near(matrixAt(exported,QStringLiteral("Sa"),2,0),1.0));
}

void TestStructuralLinearModel::unsupportedAndUnstableInputsFailClosed()
{
    RModel slack = Stage9TestModels::makeTautCable(0.0);
    slack.getProblemSetup().getStructuralCableSetup().members[0].initialTension = 0.0;
    QVERIFY_EXCEPTION_THROWN(RStructuralLinearModel::build(slack),RError);

    RModel mechanism = Stage9TestModels::makeTruss(1,0.0,0.0,false);
    QVERIFY_EXCEPTION_THROWN(RStructuralLinearModel::build(mechanism),RError);

    RModel massless = Stage9TestModels::makeTruss(1,0.0,0.0,true,false);
    QVERIFY_EXCEPTION_THROWN(RStructuralLinearModel::build(massless),RError);

    RModel timeVarying = Stage9TestModels::makeTruss();
    RConditionComponent &prescribed = timeVarying.getPoint(0).getBoundaryCondition(0)
        .getComponent(timeVarying.getPoint(0).getBoundaryCondition(0)
                          .findComponentPosition(R_VARIABLE_DISPLACEMENT_X));
    prescribed.clear();
    prescribed.add(0.0,0.0);
    prescribed.add(1.0,0.01);
    QVERIFY_EXCEPTION_THROWN(RStructuralLinearModel::build(timeVarying),RError);

    RModel mpc = Stage9TestModels::makeTruss();
    RLinearConstraint constraint;
    constraint.id = QStringLiteral("unsupported-mpc");
    constraint.terms.push_back({{1,RDofType::UX},1.0});
    mpc.getProblemSetup().getStructuralInteractionSetup().getLinearConstraints().push_back(constraint);
    QVERIFY_EXCEPTION_THROWN(RStructuralLinearModel::build(mpc),RError);

    RModel piezo = Stage9TestModels::makeTruss();
    piezo.getProblemSetup().getStructuralPiezoSetup().enabled = true;
    QVERIFY_EXCEPTION_THROWN(RStructuralLinearModel::build(piezo),RError);

    RModel thermal = Stage9TestModels::makeTruss();
    RVariable temperature(R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE);
    temperature.resize(1,thermal.getNNodes());
    for (uint node=0;node<thermal.getNNodes();++node) temperature.setValue(0,node,300.0);
    thermal.addVariable(temperature);
    QVERIFY_EXCEPTION_THROWN(RStructuralLinearModel::build(thermal),RError);
}

QTEST_APPLESS_MAIN(TestStructuralLinearModel)

#include "tst_structural_linear_model.moc"
