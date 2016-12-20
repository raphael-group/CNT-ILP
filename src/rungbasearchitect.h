#ifndef _RUNGBASEARCHITECT_H_
#define _RUNGBASEARCHITECT_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <ilcplex/ilocplex.h>
#include "copynumbertree.h"
#include "inputinstance.h"

#include "basic_types.h"
#include "basearchitect.h"

ILOSTLBEGIN

class RungBaseArchitect : public BaseArchitect
{
public:
    RungBaseArchitect(const InputInstance& inputInstance,
                      const IntMatrix &e,
                      bool rootNotFixed);

protected:
    IloIntVar4Array _a;
    IloIntVar4Array _d;
    IloIntVar4Array _A;
    IloIntVar4Array _D;

    virtual void buildVariables();
    virtual void buildConstraints();
};

#endif // _RUNGBASEARCHITECT_H_
