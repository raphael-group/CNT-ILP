#ifndef _COMPARISON_H_
#define _COMPARISON_H_

#include "copynumbertree.h"

class Comparison
{
public:
    Comparison(const CopyNumberTree& T1, const CopyNumberTree& T2);
    
    bool init();
    
    int leaf1ToLeaf2(int i) const
    {
        return _leaf1ToLeaf2[i];
    }
    
    int leaf2ToLeaf1(int j) const
    {
        return _leaf2ToLeaf1[j];
    }
    
    double robinsonFoulds() const;
    
    IntSet mapToT1(const IntSet& leafSet) const;
    
private:
    const CopyNumberTree& _T1;
    const CopyNumberTree& _T2;
    std::vector<int> _leaf1ToLeaf2;
    std::vector<int> _leaf2ToLeaf1;
};

#endif // _COMPARISON_H_
