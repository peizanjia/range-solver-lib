#ifndef REQUATIONBOOK_H
#define REQUATIONBOOK_H

#include <QList>

#include <vector>

#include <rml_dof.h>

//! Maps active nodal degrees of freedom to contiguous solver equations.
//!
//! Unlike the legacy flat node book, this class permits each node to expose a
//! different set of degrees of freedom. Equation order is deterministic:
//! nodes are visited by ID and each node's DOFs follow RDofDescriptor order.
class REquationBook
{
    private:

        uint nNodes;
        QList<RDofType> dofOrder;
        std::vector<bool> active;
        mutable std::vector<uint> equations;
        mutable uint nEquations;
        mutable bool dirty;

        uint findDofIndex(RDofType dof) const;
        uint findPosition(uint nodeId, RDofType dof) const;
        void rebuild() const;

    public:

        explicit REquationBook(uint nNodes = 0);

        //! Reset the book for the requested number of nodes; all DOFs are inactive.
        void resize(uint nNodes);

        uint getNNodes() const;

        //! Activate or deactivate a nodal DOF.
        void enable(uint nodeId, RDofType dof);
        void disable(uint nodeId, RDofType dof);

        //! Return whether a nodal DOF participates in the equation system.
        bool isEnabled(uint nodeId, RDofType dof) const;

        //! Resolve a nodal DOF to a contiguous equation number.
        bool getValue(uint nodeId, RDofType dof, uint &equation) const;

        //! Return active DOFs at a node in stable descriptor order.
        QList<RDofType> getEnabledDofs(uint nodeId) const;

        uint getNEnabled() const;
};

#endif // REQUATIONBOOK_H
