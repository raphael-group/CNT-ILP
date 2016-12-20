#ifndef _RUNGARCHITECT_H_
#define _RUNGARCHITECT_H_

#include "rungbasearchitect.h"

ILOSTLBEGIN

/// This class represent the architect that builds the ILP model for solving the
/// problem when the given input is composed of a collection of leaves and an
/// integer e representing the maximum copy-number.
class RungArchitect : public RungBaseArchitect
{
public:
    RungArchitect(const InputInstance& inputInstance,
                  const IntMatrix &e,
                  bool rootNotFixed);

protected:
    IloIntVar4Array _w;
    
    void constructTree();
    void buildVariables();
    void buildConstraints();
    void buildObjective();
};


#endif // _RUNGARCHITECT_H_
