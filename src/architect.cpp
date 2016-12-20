#include "architect.h"


Architect::Architect(const InputInstance& inputInstance,
                     const IntMatrix &e,
                     const bool rootNotFixed)
    : BaseArchitect(inputInstance, e, rootNotFixed)
    , _a()
    , _bar_a()
    , _tilde_a()
    , _d()
    , _bar_d()
    , _tilde_d()
    , _w()
{

}

void Architect::constructTree()
{
    // construct tree topology
    BaseArchitect::constructTree();
    
    // construct events
    const CopyNumberTree::Digraph& T = _T.T();
    for (CopyNumberTree::ArcIt a_ij(T); a_ij != lemon::INVALID; ++a_ij)
    {
        CopyNumberTree::Node v_i = T.source(a_ij);
        CopyNumberTree::Node v_j = T.target(a_ij);
        const int i = _T.index(v_i);
        const int j = _T.index(v_j);
        const int remapped_j = remap_j(i, j);
        
        for (int chr = 0; chr < _num_chr; ++chr)
        {
            for (int s = 0; s < _n[chr]; ++s)
            {
                for (int t = s; t < _n[chr]; ++t)
                {
                    const int remapped_t = remap_t(chr, s, t);
                    const int d_chr_i_j_s_t = _cplex.getIntValue(_d[chr][i][remapped_j][s][remapped_t]);
                    if (d_chr_i_j_s_t > 0)
                    {
                        _T.addEvent(chr, i, j, CopyNumberTree::Event(chr, s, t, -d_chr_i_j_s_t));
                    }
                }
            }
            
            for (int s = 0; s < _n[chr]; ++s)
            {
                for (int t = s; t < _n[chr]; ++t)
                {
                    const int remapped_t = remap_t(chr, s, t);
                    const int a_chr_i_j_s_t = _cplex.getIntValue(_a[chr][i][remapped_j][s][remapped_t]);
                    if (a_chr_i_j_s_t > 0)
                    {
                        _T.addEvent(chr, i, j, CopyNumberTree::Event(chr, s, t, a_chr_i_j_s_t));
                    }
                }
            }
        }
    }
}

void Architect::buildVariables()
{
    BaseArchitect::buildVariables();
    
    char buf[1024];
    
    _a = IloIntVar5Array(_env, _num_chr);
    _d = IloIntVar5Array(_env, _num_chr);
    _w = IloIntVar5Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _a[chr] = IloIntVar4Array(_env, _k - 1);
        _d[chr] = IloIntVar4Array(_env, _k - 1);
        _w[chr] = IloIntVar4Array(_env, _k - 1);
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            _a[chr][i] = IloIntVar3Array(_env, _num_vertices - (i + 1));
            _d[chr][i] = IloIntVar3Array(_env, _num_vertices - (i + 1));
            _w[chr][i] = IloIntVar3Array(_env, _num_vertices - (i + 1));
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                _a[chr][i][remapped_j] = IloIntVarMatrix(_env, _n[chr]);
                _d[chr][i][remapped_j] = IloIntVarMatrix(_env, _n[chr]);
                _w[chr][i][remapped_j] = IloIntVarMatrix(_env, _n[chr]);
                for(unsigned int s = 0; s < _n[chr]; ++s)
                {
                    _a[chr][i][remapped_j][s] = IloIntVarArray(_env, _n[chr] - s);
                    _d[chr][i][remapped_j][s] = IloIntVarArray(_env, _n[chr] - s);
                    _w[chr][i][remapped_j][s] = IloIntVarArray(_env, _n[chr] - s);
                    for(unsigned int t = s; t < _n[chr]; ++t)
                    {
                        const int remapped_t = remap_t(chr, s, t);
                        const unsigned int min_e = min(_e[chr][s], _e[chr][t]);
                        
                        snprintf(buf, 1024, "a_%d_%d_%d_%d_%d", chr, i, j, s, t);
                        _a[chr][i][remapped_j][s][remapped_t] = IloIntVar(_env, 0, min_e, buf);
                        
                        snprintf(buf, 1024, "d_%d_%d_%d_%d_%d", chr, i, j, s, t);
                        _d[chr][i][remapped_j][s][remapped_t] = IloIntVar(_env, 0, min_e, buf);
                        
                        snprintf(buf, 1024, "w_%d_%d_%d_%d_%d", chr, i, j, s, t);
                        _w[chr][i][remapped_j][s][remapped_t] = IloIntVar(_env, 0, min_e + min_e, buf);
                    }
                }
            }
        }
    }
  
    _bar_a = IloIntVar4Array(_env, _num_chr);
    _tilde_a = IloIntVar4Array(_env, _num_chr);
    _bar_d = IloIntVar4Array(_env, _num_chr);
    _tilde_d = IloIntVar4Array(_env, _num_chr);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        _bar_a[chr] = IloIntVar3Array(_env, _k - 1);
        _tilde_a[chr] = IloIntVar3Array(_env, _k - 1);
        _bar_d[chr] = IloIntVar3Array(_env, _k - 1);
        _tilde_d[chr] = IloIntVar3Array(_env, _k - 1);
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            _bar_a[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            _tilde_a[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            _bar_d[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            _tilde_d[chr][i] = IloIntVarMatrix(_env, _num_vertices - (i + 1));
            for(unsigned j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                _bar_a[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                _tilde_a[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                _bar_d[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                _tilde_d[chr][i][remapped_j] = IloIntVarArray(_env, _n[chr]);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    const unsigned int max_e = *std::max_element(_e[chr].begin(), _e[chr].end());
                    
                    snprintf(buf, 1024, "bar_a_%d_%d_%d_%d", chr, i, j, l);
                    _bar_a[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr][l], buf);
                    
                    snprintf(buf, 1024, "tilde_a_%d_%d_%d_%d", chr, i, j, l);
                    _tilde_a[chr][i][remapped_j][l] = IloIntVar(_env, 0, max_e, buf);
                    
                    snprintf(buf, 1024, "bar_d_%d_%d_%d_%d", chr, i, j, l);
                    _bar_d[chr][i][remapped_j][l] = IloIntVar(_env, 0, _e[chr][l], buf);
                    
                    snprintf(buf, 1024, "tilde_d_%d_%d_%d_%d", chr, i, j, l);
                    _tilde_d[chr][i][remapped_j][l] = IloIntVar(_env, 0, max_e, buf);
                }
            }
        }
    }    
}

void Architect::buildConstraints()
{
    BaseArchitect::buildConstraints();
    
    char buf[1024];

    /**
       1. y_{j,\ell} = y_{i,\ell} - \bar{d}_{i,j,\ell} + \bar{a}_{i,j,\ell} for 1 <= \ell <= n and (v_i,v_j) \in E(G)
       2. \sum_{s <= \ell <= t} a_{i,j,s,t} = \bar{a}_{i,j,\ell} + \tilde{a}_{i,j,\ell} and 1 <= \ell <= n, (v_i,v_j) \in E(G)
       3. \sum_{s <= \ell <= t} d_{i,j,s,t} = \bar{d}_{i,j,\ell} + \tilde{d}_{i,j,\ell} and 1 <= \ell <= n, (v_i,v_j) \in E(G)
     **/
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned int j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int l = 0; l < _n[chr]; ++l)
                {
                    IloConstraint cons(_y[chr][j][l]
                                       - _y[chr][i][l]
                                       + _bar_d[chr][i][remapped_j][l]
                                       - _bar_a[chr][i][remapped_j][l] == 0);
                    snprintf(buf, 1024, "apply_real_events_%d_%d_%d_%d", chr, i, j, l);
                    cons.setName(buf);
                    _model.add(cons);
                    
                    IloExpr sum_amp(_env);
                    IloExpr sum_del(_env);
                    for(unsigned int s = 0; s <= l; ++s)
                    {
                        for(unsigned int t = l; t < _n[chr]; ++t)
                        {
                            const int remapped_t = remap_t(chr, s, t);
                            sum_amp += _a[chr][i][remapped_j][s][remapped_t];
                            sum_del += _d[chr][i][remapped_j][s][remapped_t];
                        }
                    }

                    IloConstraint cons2(sum_amp - _bar_a[chr][i][remapped_j][l] - _tilde_a[chr][i][remapped_j][l] == 0);
                    snprintf(buf, 1024, "amp_%d_%d_%d_%d", chr, i, j, l);
                    cons2.setName(buf);
                    _model.add(cons2);
                    
                    IloConstraint cons3(sum_del - _bar_d[chr][i][remapped_j][l] - _tilde_d[chr][i][remapped_j][l] == 0);
                    snprintf(buf, 1024, "del_%d_%d_%d_%d", chr, i, j, l);
                    cons3.setName(buf);
                    _model.add(cons3);
                }
            }
        }
    }

    /**
       1. \bar{a}_{i,j,s} <= e \cdot \bar{y}_{j, s} for 1 <= s <= n and (v_i,v_j) \in E(G)
       2. \tilde{a}_{i,j,s} <= e \cdot (1 - \bar{y}_{j, s}) for 1 <= s <= n and (v_i,v_j) \in E(G)
       3. \bar{d}_{i,j,s} <= e \cdot (\bar{y}_{i, s} + \bar{y}_{j, s}) for 1 <= s <= n and (v_i,v_j) \in E(G)
       4. \tilde{d}_{i,j,s} <= e \cdot (2 - \bar{y}_{i,s} - \bar{y}_{j,s}) for 1 <= s <= n and (v_i,v_j) \in E(G)
       5. \bar{d}_{i, j, s} <= y_{i,s} - \bar{y}_{j,s} for 1 <= s <= n and (v_i,v_j) \in E(G)
     **/
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned int j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int s = 0; s < _n[chr]; ++s)
                {
                    const int max_e = *max_element(_e[chr].begin(), _e[chr].end());
                    const int e_chr_s = _e[chr][s];

                    IloConstraint cons(_bar_a[chr][i][remapped_j][s] - e_chr_s * _bar_y[chr][j][s] <= 0);
                    snprintf(buf, 1024, "real_amp_%d_%d_%d_%d", chr, i, j, s);
                    cons.setName(buf);
                    _model.add(cons);
                    
                    IloConstraint cons2(_tilde_a[chr][i][remapped_j][s] - max_e + max_e * _bar_y[chr][j][s] <= 0);
                    snprintf(buf, 1024, "fake_amp_%d_%d_%d_%d", chr, i, j, s);
                    cons2.setName(buf);
                    _model.add(cons2);
                    
                    IloConstraint cons3(_bar_d[chr][i][remapped_j][s] - e_chr_s * _bar_y[chr][i][s] - e_chr_s * _bar_y[chr][j][s] <= 0);
                    snprintf(buf, 1024, "real_del_%d_%d_%d_%d", chr, i, j, s);
                    cons3.setName(buf);
                    _model.add(cons3);
                    
                    IloConstraint cons4(_tilde_d[chr][i][remapped_j][s] - max_e * 2 + max_e * _bar_y[chr][i][s] + max_e * _bar_y[chr][j][s] <= 0);
                    snprintf(buf, 1024, "fake_del_%d_%d_%d_%d", chr, i, j, s);
                    cons4.setName(buf);
                    _model.add(cons4);
                    
                    IloConstraint cons5(_bar_d[chr][i][remapped_j][s] - _y[chr][i][s] + _bar_y[chr][j][s] - 1 + _x[i][remapped_j] <= 0);
                    snprintf(buf, 1024, "feasibility_%d_%d_%d_%d", chr, i, j, s);
                    cons5.setName(buf);
                    _model.add(cons5);
                }
            }
        }
    }

    /**
       w_{i,j,s,t} >= a_{i,j,s,t} + d_{i,j,s,t} - (1 - x_{i,j}) \cdot 2e for 1 <= s <= t <= n and (v_i,v_j) \in E(G)
     **/
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned int j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int s = 0; s < _n[chr]; ++s)
                {
                    for(unsigned int t = s; t < _n[chr]; ++t)
                    {
                        const int remapped_t = remap_t(chr, s, t);
                        const int min_e = min(_e[chr][s], _e[chr][t]);
                        IloConstraint cons(_w[chr][i][remapped_j][s][remapped_t]
                                           - _a[chr][i][remapped_j][s][remapped_t]
                                           - _d[chr][i][remapped_j][s][remapped_t]
                                           + 2*min_e - (2*min_e) * _x[i][remapped_j] >= 0);
                        snprintf(buf, 1024, "cost_%d_%d_%d_%d_%d", chr, i, j, s, t);
                        cons.setName(buf);
                        _model.add(cons);
                    }
                }
            }
        }
    }
}

void Architect::buildObjective()
{
    _obj = IloExpr(_env);
    for(unsigned int chr = 0; chr < _num_chr; ++chr)
    {
        for(unsigned int i = 0; i < (_k - 1); ++i)
        {
            for(unsigned int j = (i + 1); j < _num_vertices; ++j)
            {
                const int remapped_j = remap_j(i, j);
                for(unsigned int s = 0; s < _n[chr]; ++s)
                {
                    for(unsigned int t = s; t < _n[chr]; ++t)
                    {
                        const int remapped_t = remap_t(chr, s, t);
                        _obj += _w[chr][i][remapped_j][s][remapped_t];
                    }
                }
            }
        }
    }
    _model.add(IloMinimize(_env, _obj));
}
