#ifndef _BASEARCHITECT_H_
#define _BASEARCHITECT_H_

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

ILOSTLBEGIN

class BaseArchitect
{
public:
    BaseArchitect(const InputInstance& inputInstance,
                  const IntMatrix &e,
                  bool rootNotFixed);
    virtual ~BaseArchitect()
    {
    }
    
    void init();
    bool solve(const int timeLimit, const int memoryLimit);
    
    double getOpt()
    {
        return _cplex.getObjValue();
    }
    
    const CopyNumberTree& getTree()
    {
        return _T;
    }
    
    void exportModel(const std::string& filename) const;
    
protected:
    typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
    typedef IloArray<IloBoolVarMatrix> IloBoolVar3Array;
    typedef IloArray<IloBoolVar3Array> IloBoolVar4Array;
    typedef IloArray<IloBoolVar4Array> IloBoolVar5Array;
    
    typedef IloArray<IloNumVarArray> IloNumVarMatrix;
    typedef IloArray<IloNumVarMatrix> IloNumVar3Array;
    typedef IloArray<IloNumVar3Array> IloNumVar4Array;
    typedef IloArray<IloNumVar4Array> IloNumVar5Array;
    
    typedef IloArray<IloIntVarArray> IloIntVarMatrix;
    typedef IloArray<IloIntVarMatrix> IloIntVar3Array;
    typedef IloArray<IloIntVar3Array> IloIntVar4Array;
    typedef IloArray<IloIntVar4Array> IloIntVar5Array;
    
    int remap_j(int i, int j) const
    {
        assert(0 <= i && i < _k - 1);
        return j - (i + 1);
    }
    
    virtual void constructTree();
    
protected:
    const Int3Array& _C;
    const IntMatrix& _e;
    const bool _rootNotFixed;
    const unsigned int _num_chr;
    const unsigned int _k;
    const unsigned int _num_vertices;
    const IntArray& _n;

    IloEnv _env;
    IloModel _model;
    IloCplex _cplex;
    CopyNumberTree _T;
    
    /// _x[i][j] == 1 if and only if there is an edge (i,j)
    IloBoolVarMatrix _x;
    IloIntVar3Array _y;
    IloBoolVar3Array _bar_y;
    IntMatrix _num_z;
    IloBoolVar4Array _z;
    
    IloExpr _obj;
    
    virtual void buildVariables();
    virtual void buildConstraints();
    virtual void buildObjective() = 0;
};

#endif // _BASEARCHITECT_H_
