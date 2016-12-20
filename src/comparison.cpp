#include "comparison.h"

Comparison::Comparison(const CopyNumberTree& T1, const CopyNumberTree& T2)
    : _T1(T1)
    , _T2(T2)
{
}

bool Comparison::init()
{
    // check leaves
    if (_T1.k() != _T2.k())
    {
        return false;
    }
    
    // check #chromosomes
    if (_T1.numChr() != _T2.numChr())
    {
        return false;
    }
    
    // check #segments per chromosome
    const IntArray& n1 = _T1.n();
    const IntArray& n2 = _T2.n();
    for (int chr = 0; chr < _T1.numChr(); ++chr)
    {
        if (n1[chr] != n2[chr])
        {
            return false;
        }
    }
    
    // construct leaf to leaf mapping
    const int k = _T1.k();
    
    _leaf1ToLeaf2 = std::vector<int>(2*k - 1, -1);
    _leaf2ToLeaf1 = std::vector<int>(2*k - 1, -1);
    
    for (int i = k - 1; i < 2*k - 1; ++i)
    {
        for (int j = k - 1; j < 2*k - 1; ++j)
        {
            if (_T1.profile(i) == _T2.profile(j))
            {
                _leaf1ToLeaf2[i] = j;
                _leaf2ToLeaf1[j] = i;
            }
        }
    }
    
    // do we have a mapping?
    for (int i = k-1; i < 2*k-1; ++i)
    {
        if (_leaf1ToLeaf2[i] == -1)
        {
            return false;
        }
        if (_leaf2ToLeaf1[i] == -1)
        {
            return false;
        }
    }
    
    return true;
}

double Comparison::robinsonFoulds() const
{
    std::set<IntSetPair> T1_splits;
    std::set<IntSetPair> T2_splits;
    
    const int k = _T1.k();
    for (int i = 0; i < k - 1; ++i)
    {
        T1_splits.insert(_T1.splits(i));

        IntSetPair split = _T2.splits(i);
        split.first = mapToT1(split.first);
        split.second = mapToT1(split.second);
        T2_splits.insert(split);
    }
    
    std::set<IntSetPair> result;
    std::set_symmetric_difference(T1_splits.begin(), T1_splits.end(),
                                  T2_splits.begin(), T2_splits.end(),
                                  std::inserter(result, result.begin()));
    
    return result.size() / (2.0*(k-1));
}

IntSet Comparison::mapToT1(const IntSet& leafSet) const
{
    IntSet result;
    
    for (int i : leafSet)
    {
        assert(_leaf2ToLeaf1[i] != -1);
        result.insert(_leaf2ToLeaf1[i]);
    }
    
    return result;
}
