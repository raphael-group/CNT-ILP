#ifndef _ARCHITECT_H_
#define _ARCHITECT_H_

#include "basearchitect.h"

ILOSTLBEGIN

/// This class represent the architect that builds the ILP model for solving the
/// problem when the given input is composed of a collection of leaves and an
/// integer e representing the maximum copy-number.
class Architect : public BaseArchitect
{
public:
    Architect(const InputInstance& inputInstance,
              const IntMatrix &e,
              bool rootNotFixed);
  
protected:
    int remap_t(int chr, int s, int t) const
    {
        assert(0 <= s && s < _n[chr] && s <= t && t < _n[chr]);
        return t - s;
    }
    
    void constructTree();

protected:
    IloIntVar5Array _a;
    IloIntVar4Array _bar_a;
    IloIntVar4Array _tilde_a;
    IloIntVar5Array _d;
    IloIntVar4Array _bar_d;
    IloIntVar4Array _tilde_d;
    IloIntVar5Array _w;

    void buildVariables();
    void buildConstraints();
    void buildObjective();
};

#endif // _ARCHITECT_H_
