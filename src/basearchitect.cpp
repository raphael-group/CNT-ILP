#include "basearchitect.h"
#include <lemon/time_measure.h>

BaseArchitect::BaseArchitect(const InputInstance& inputInstance,
                             const IntMatrix &e,
                             bool rootNotFixed)
    : _C(inputInstance.C())
    , _e(e)
    , _rootNotFixed(rootNotFixed)
    , _num_chr(inputInstance.numChr())
    , _k(inputInstance.k())
    , _num_vertices(2*_k - 1)
    , _n(inputInstance.n())
    , _env()
    , _model(_env)
    , _cplex(_model)
    , _T(inputInstance.k(), inputInstance.numChr(), inputInstance.n())
    , _x()
    , _y()
    , _bar_y()
    , _num_z()
    , _z()
    , _obj()
{
}

void BaseArchitect::init()
{
    buildVariables();
    buildConstraints();
    buildObjective();
}

bool BaseArchitect::solve(const int timeLimit, const int memoryLimit)
{
    _cplex.setOut(std::cerr);
    _cplex.setWarning(std::cerr);
    _cplex.setError(std::cerr);
    
    // shut up cplex
    bool verbose = true;
    if (!verbose) {
        _cplex.setOut(_env.getNullStream());
        _cplex.setWarning(_env.getNullStream());
        _cplex.setError(_env.getNullStream());
    }
    
    if (timeLimit > 0)
    {
        _cplex.setParam(IloCplex::TiLim, timeLimit);
    }
    
    if (memoryLimit > 0)
    {
        _cplex.setParam(IloCplex::TreLim, memoryLimit);
    }
    
    lemon::Timer timer;
    bool res = _cplex.solve();
    if (res)
    {
        _cplex.out() << std::endl;
        _cplex.out() << "Solution status = " << _cplex.getStatus() << std::endl;
        _cplex.out() << "Solution LB = " << _cplex.getBestObjValue() << std::endl;
        _cplex.out() << "Solution UB = " << _cplex.getObjValue() << std::endl;
        _cplex.out() << "Runtime = " << timer.realTime() << " seconds" << std::endl;
        
        constructTree();
    }
    else
    {
        _cplex.out() << "FAILED TO SOLVE:  " << _cplex.getStatus() << std::endl;
    }
        
    return res;
}

void BaseArchitect::exportModel(const std::string& filename) const
{
    _cplex.exportModel(filename.c_str());
}

void BaseArchitect::constructTree()
{
    // add and label arcs
    for(unsigned int i = 0; i < (_k - 1); ++i)
    {
        for(unsigned j = (i + 1); j < _num_vertices; ++j)
        {
            const int remapped_j = remap_j(i,j);
            if(_cplex.getIntValue(_x[i][remapped_j]))
            {
                _T.addArc(i, j);
            }
        }
    }
    
    // construct profiles
    for(unsigned int i = 0; i < _num_vertices; ++i)
    {
        CopyNumberTree::ProfileVector profile(_num_chr);
        for (int chr = 0; chr < _num_chr; ++chr)
        {
            profile[chr] = CopyNumberTree::Profile(_n[chr], 0);
            for (int s = 0; s < _n[chr]; ++s)
            {
                profile[chr][s] = _cplex.getIntValue(_y[chr][i][s]);
            }
        }
        _T.setProfile(i, profile);
    }
}

void BaseArchitect::buildVariables()
{
    char buf[1024];
    
    _x = IloBoolVarMatrix(_env, _k - 1);
    for(unsigned int i = 0; i < (_k - 1); ++i)
    {
        _x[i] = IloBoolVarArray(_env, _num_vertices - (i + 1));
        for(unsigned j = (i + 1); j < _num_vertices; ++j)
        {
            snprintf(buf, 1024, "x_%d_%d", i, j);
            _x[i][remap_j(i, j)] = IloBoolVar(_env, buf);
        }
    }
    
    _y = IloIntVar3Array(_env, _num_chr);
    _bar_y = IloBoolVar3Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _y[chr] = IloIntVarMatrix(_env, _num_vertices);
        _bar_y[chr] = IloBoolVarMatrix(_env, _num_vertices);
        for(unsigned int i = 0; i < _num_vertices; ++i)
        {
            _y[chr][i] = IloIntVarArray(_env, _n[chr]);
            _bar_y[chr][i] = IloBoolVarArray(_env, _n[chr]);
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                snprintf(buf, 1024, "y_%d_%d_%d", chr, i, s);
                _y[chr][i][s] = IloIntVar(_env, 0, _e[chr][s], buf);
                
                snprintf(buf, 1024, "bar_y_%d_%d_%d", chr, i, s);
                _bar_y[chr][i][s] = IloBoolVar(_env, buf);
            }
        }
    }
    
    _num_z = IntMatrix(_num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _num_z[chr] = IntArray(_n[chr]);
        for(unsigned int s = 0; s < _n[chr]; ++s)
        {
            _num_z[chr][s] = floor(log(_e[chr][s]) / log(2.0)) + 1;
        }
    }
    
    _z = IloBoolVar4Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _z[chr] = IloBoolVar3Array(_env, _num_vertices);
        for(unsigned int i = 0; i < _num_vertices; ++i)
        {
            _z[chr][i] = IloBoolVarMatrix(_env, _n[chr]);
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                _z[chr][i][s] = IloBoolVarArray(_env, _num_z[chr][s]);
                for(unsigned int q = 0; q < _num_z[chr][s]; ++q)
                {
                    snprintf(buf, 1024, "z_%d_%d_%d_%d", chr, i, s, q);
                    _z[chr][i][s][q] = IloBoolVar(_env, buf);
                }
            }
        }
    }
}

void BaseArchitect::buildConstraints()
{
    char buf[1024];
    /**
     \sum_{i \in \delta^-(j)} x_{i,j} = 1 for 1 < j <= 2k-1
     **/
    for(unsigned int j = 1; j < _num_vertices; ++j)
    {
        IloExpr sum_x_minus(_env);
        unsigned int lim = std::min(j, _k - 1);
        
        for(unsigned i = 0; i < lim; ++i)
        {
            sum_x_minus += _x[i][remap_j(i, j)];
        }
        
        IloConstraint cons(sum_x_minus == 1);
        snprintf(buf, 1024, "in_deg_%d", j);
        cons.setName(buf);
        _model.add(cons);
    }
    
    /**
     \sum_{j \in \delta^+(i)} x_{i,j} = 2 for 1 <= i < k
     **/
    for(unsigned int i = 0; i < (_k - 1); ++i)
    {
        IloExpr sum_x_plus(_env);
        
        for(unsigned int j = (i + 1); j < _num_vertices; ++j)
        {
            sum_x_plus += _x[i][remap_j(i,j)];
        }
        
        IloConstraint cons(sum_x_plus == 2);
        snprintf(buf, 1024, "out_deg_%d", i);
        cons.setName(buf);
        _model.add(cons);
    }
    
    /**
     y_{1,s} = 2 for 1 <= s <= n
     **/
    if (!_rootNotFixed)
    {
        for(unsigned int chr = 0; chr < _num_chr; ++chr)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                IloConstraint cons(_y[chr][0][s] == 2);
                snprintf(buf, 1024, "root_%d", s);
                cons.setName(buf);
                _model.add(cons);
            }
        }
    }
    
    /**
     y_{i,s} = c_{i-k+2,s} for  k <= i <= 2k-1 and 1 <= s <= n
     **/
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = (_k - 1); i < _num_vertices; ++i)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                IloConstraint cons(_y[chr][i][s] == _C[chr][i - _k + 1][s]);
                snprintf(buf, 1024, "leaf_%d_%d", i, s);
                cons.setName(buf);
                _model.add(cons);
            }
        }
    }
    
    /**
     1. y_{i,s} = \sum_{j=0}^{\lceil \log_2(e) \rceil} 2^j \cdot z_{i,s,j} for 1 <= i <= 2k-1 and 1 <= s <= n
     2. \bar{y}_{i,s}  <= \sum_{j=0}^{\lceil \log_2(e) \rceil} z_{i,s,j} for 1 <= i <= 2k-1 and 1 <= s <= n
     3. \bar{y}_{i,s}  \ge z_{i,s,j} for  1 <= i <= 2k-1 and 1 <= s <= n and 0 <= j <= (\lceil \log_2(e) \rceil)
     **/
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < _num_vertices; ++i)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                IloExpr sum_z_2(_env);
                IloExpr sum_z(_env);
                
                for(unsigned int q = 0; q < _num_z[chr][s]; ++q)
                {
                    sum_z_2 += (1 << q) * _z[chr][i][s][q];
                    sum_z += _z[chr][i][s][q];
                    
                    IloConstraint cons(_bar_y[chr][i][s] - _z[chr][i][s][q] >= 0);
                    snprintf(buf, 1024, "nonzero_lb_%d_%d_%d_%d", chr, i, s, q);
                    cons.setName(buf);
                    _model.add(cons);
                }
                
                IloConstraint cons2(_y[chr][i][s] - sum_z_2 == 0);
                snprintf(buf, 1024, "binary_%d_%d_%d", chr, i, s);
                cons2.setName(buf);
                _model.add(cons2);
                
                IloConstraint cons3(_bar_y[chr][i][s] - sum_z <= 0);
                snprintf(buf, 1024, "nonzero_ub_%d_%d_%d", chr, i, s);
                cons3.setName(buf);
                _model.add(cons3);
            }
        }
    }
 }
