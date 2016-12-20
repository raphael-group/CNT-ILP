#include "rungbasearchitect.h"

RungBaseArchitect::RungBaseArchitect(const InputInstance& inputInstance,
                                     const IntMatrix &e,
                                     bool rootNotFixed)
    : BaseArchitect(inputInstance, e, rootNotFixed)
    , _a()
    , _d()
    , _A()
    , _D()
{
}

void RungBaseArchitect::buildVariables()
{
    BaseArchitect::buildVariables();
    
    _a = IloIntVar4Array(_env, _num_chr);
    _d = IloIntVar4Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        const int e = *max_element(_e[chr].begin(), _e[chr].end());
        _a[chr] = IloIntVar3Array(_env, _k - 1);
        _d[chr] = IloIntVar3Array(_env, _k - 1);
        for(unsigned int i = 0; i < _k - 1; ++i)
        {
            _a[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            _d[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            for(unsigned int j = i + 1; j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                _a[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                _d[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                for(unsigned int s = 0; s < _n[chr]; ++s)
                {
                    _a[chr][i][remapped_j][s] = IloIntVar(_env, 0, e);
                    _d[chr][i][remapped_j][s] = IloIntVar(_env, 0, e);
                }
            }
        }
    }
    
    _A = IloIntVar4Array(_env, _num_chr);
    _D = IloIntVar4Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _A[chr] = IloIntVar3Array(_env, _k - 1);
        _D[chr] = IloIntVar3Array(_env, _k - 1);
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            _A[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            _D[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                _A[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                _D[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    _A[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr][l]);
                    _D[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr][l]);
                }
            }
        }
    }
}

void RungBaseArchitect::buildConstraints()
{
    BaseArchitect::buildConstraints();
    
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    if(l == 0)
                    {
                        _model.add(_A[chr][i][remapped_j][l] - _a[chr][i][remapped_j][l] >= 0);
                        _model.add(_D[chr][i][remapped_j][l] - _d[chr][i][remapped_j][l] >= 0);
                    }
                    else
                    {
                        _model.add(_A[chr][i][remapped_j][l] - _a[chr][i][remapped_j][l] + _a[chr][i][remapped_j][l-1] >= 0);
                        _model.add(_D[chr][i][remapped_j][l] - _d[chr][i][remapped_j][l] + _d[chr][i][remapped_j][l-1] >= 0);
                    }
                }
            }
        }
    }
}
