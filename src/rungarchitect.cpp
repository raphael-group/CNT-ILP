#include "rungarchitect.h"

RungArchitect::RungArchitect(const InputInstance& inputInstance,
                                     const IntMatrix &e,
                                     bool rootNotFixed)
    : RungBaseArchitect(inputInstance, e, rootNotFixed)
    , _w()
{
}

void RungArchitect::constructTree()
{
    // construct tree topology
    BaseArchitect::constructTree();
    
    //collect the number of amplifications and deletions per segment
    Int4Array amplifications(_num_chr);
    Int4Array deletions(_num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        amplifications[chr] = Int3Array(_k - 1);
        deletions[chr] = Int3Array(_k - 1);
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            amplifications[chr][i] = IntMatrix(_num_vertices - (i + 1));
            deletions[chr][i] = IntMatrix(_num_vertices - (i + 1));
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                amplifications[chr][i][remapped_j] = IntArray(_n[chr]);
                deletions[chr][i][remapped_j] = IntArray(_n[chr]);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    amplifications[chr][i][remapped_j][l] = _cplex.getIntValue(_a[chr][i][remapped_j][l]);
                    deletions[chr][i][remapped_j][l] = _cplex.getIntValue(_d[chr][i][remapped_j][l]);
                }
            }
        }
    }
    
    // construct events
    const CopyNumberTree::Digraph& T = _T.T();
    for (CopyNumberTree::ArcIt a_ij(T); a_ij != lemon::INVALID; ++a_ij)
    {
        CopyNumberTree::Node v_i = T.source(a_ij);
        CopyNumberTree::Node v_j = T.target(a_ij);
        const int i = _T.index(v_i);
        const int j = _T.index(v_j);
        const int remapped_j = remap_j(i, j);
        unsigned int counter = 0;
        
        for (unsigned int chr = 0; chr < _num_chr; ++chr)
        {
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                while(deletions[chr][i][remapped_j][s] != 0)
                {
                    int t = s;
                    while (t+1 < _n[chr] && deletions[chr][i][remapped_j][t+1] != 0) ++t;
                    
                    _T.addEvent(chr, i, j, CopyNumberTree::Event(chr, s, t, -1));
                    ++counter;
                    
                    for (int l = s; l <= t; ++l)
                        --deletions[chr][i][remapped_j][l];
                }
            }
            
            for(unsigned int s = 0; s < _n[chr]; ++s)
            {
                while(amplifications[chr][i][remapped_j][s] != 0)
                {
                    int t = s;
                    while (t+1 < _n[chr] && amplifications[chr][i][remapped_j][t+1] != 0) ++t;
                    
                    _T.addEvent(chr, i, j, CopyNumberTree::Event(chr, s, t, 1));
                    ++counter;
                    
                    for (int l = s; l <= t; ++l)
                        --amplifications[chr][i][remapped_j][l];
                }
            }
            
            unsigned int w = 0;
            for(unsigned int l = 0; l < _n[chr]; ++l)
            {
                w += _cplex.getIntValue(_w[chr][i][remapped_j][l]);
            }
            
            assert(counter == w);
        }
    }
}


void RungArchitect::buildVariables()
{
    RungBaseArchitect::buildVariables();
    
    _w = IloIntVar4Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _w[chr] = IloIntVar3Array(_env, _k - 1);
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            _w[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                _w[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    _w[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr][l] + _e[chr][l]);
                }
            }
        }
    }
}

void RungArchitect::buildConstraints()
{
    RungBaseArchitect::buildConstraints();
    
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        const int e = *max_element(_e[chr].begin(), _e[chr].end());
        for(unsigned int i = 0; i < _k - 1; ++i)
        {
            for(unsigned int j = i + 1; j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int s = 0; s < _n[chr]; ++s)
                {
                    _model.add(_bar_y[chr][i][s] - _bar_y[chr][j][s] + 1 - _x[i][remapped_j] >= 0);
                    _model.add(_y[chr][i][s]
                               - _d[chr][i][remapped_j][s]
                               - e
                               + e * _bar_y[chr][i][s]
                               - e * _bar_y[chr][j][s] <= 0);
                    _model.add(_y[chr][j][s]
                               - _y[chr][i][s]
                               + _d[chr][i][remapped_j][s]
                               - _a[chr][i][remapped_j][s]
                               - 4 * e
                               + (2 * e) * _bar_y[chr][i][s]
                               + (2 * e) * _bar_y[chr][j][s] <= 0);
                    _model.add(_y[chr][j][s]
                               - _y[chr][i][s]
                               + _d[chr][i][remapped_j][s]
                               - _a[chr][i][remapped_j][s]
                               + 4 * e
                               - (2 * e) * _bar_y[chr][i][s]
                               - (2 * e) * _bar_y[chr][j][s] >= 0);
                    _model.add(_d[chr][i][remapped_j][s]
                               - _y[chr][i][s]
                               + 1
                               - (e + 1) * 2
                               + (e + 1) * _bar_y[chr][i][s]
                               + (e + 1) * _bar_y[chr][j][s] <= 0);
                    
                    _model.add(_w[chr][i][remapped_j][s]
                               - _A[chr][i][remapped_j][s]
                               - _D[chr][i][remapped_j][s]
                               + (2 * e)
                               - (2 * e) * _x[i][remapped_j] >= 0);
                }
            }
        }
    }
}

void RungArchitect::buildObjective()
{
    _obj = IloExpr(_env);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    _obj += _w[chr][i][remapped_j][l];
                }
            }
        }
    }
    _model.add(IloMinimize(_env, _obj));
}
