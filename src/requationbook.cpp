#include <rbl_error.h>
#include <rbl_utils.h>

#include "requationbook.h"

REquationBook::REquationBook(uint nNodes)
    : nNodes(0)
    , nEquations(0)
    , dirty(false)
{
    const QList<RDofDescriptor> descriptors = RDofDescriptor::getAll();
    this->dofOrder.reserve(descriptors.size());
    for (const RDofDescriptor &descriptor : descriptors)
    {
        this->dofOrder.append(descriptor.getType());
    }
    this->resize(nNodes);
}

void REquationBook::resize(uint nNodes)
{
    this->nNodes = nNodes;
    const size_t size = size_t(nNodes)*size_t(this->dofOrder.size());
    this->active.assign(size,false);
    this->equations.assign(size,RConstants::eod);
    this->nEquations = 0;
    this->dirty = false;
}

uint REquationBook::getNNodes() const
{
    return this->nNodes;
}

void REquationBook::enable(uint nodeId, RDofType dof)
{
    const uint position = this->findPosition(nodeId,dof);
    if (!this->active[position])
    {
        this->active[position] = true;
        this->dirty = true;
    }
}

void REquationBook::disable(uint nodeId, RDofType dof)
{
    const uint position = this->findPosition(nodeId,dof);
    if (this->active[position])
    {
        this->active[position] = false;
        this->dirty = true;
    }
}

bool REquationBook::isEnabled(uint nodeId, RDofType dof) const
{
    return this->active[this->findPosition(nodeId,dof)];
}

bool REquationBook::getValue(uint nodeId, RDofType dof, uint &equation) const
{
    const uint position = this->findPosition(nodeId,dof);
    if (!this->active[position])
    {
        equation = RConstants::eod;
        return false;
    }
    this->rebuild();
    equation = this->equations[position];
    return true;
}

QList<RDofType> REquationBook::getEnabledDofs(uint nodeId) const
{
    R_ERROR_ASSERT(nodeId < this->nNodes);
    QList<RDofType> dofs;
    for (RDofType dof : this->dofOrder)
    {
        if (this->isEnabled(nodeId,dof))
        {
            dofs.append(dof);
        }
    }
    return dofs;
}

uint REquationBook::getNEnabled() const
{
    this->rebuild();
    return this->nEquations;
}

uint REquationBook::findDofIndex(RDofType dof) const
{
    const qsizetype index = this->dofOrder.indexOf(dof);
    R_ERROR_ASSERT(index >= 0);
    return uint(index);
}

uint REquationBook::findPosition(uint nodeId, RDofType dof) const
{
    R_ERROR_ASSERT(nodeId < this->nNodes);
    return nodeId*uint(this->dofOrder.size()) + this->findDofIndex(dof);
}

void REquationBook::rebuild() const
{
    if (!this->dirty)
    {
        return;
    }
    uint nextEquation = 0;
    for (size_t i=0;i<this->active.size();i++)
    {
        this->equations[i] = this->active[i] ? nextEquation++ : RConstants::eod;
    }
    this->nEquations = nextEquation;
    this->dirty = false;
}
