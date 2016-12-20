#ifndef _INPUTINSTANCE_H
#define _INPUTINSTANCE_H

#include "basic_types.h"
#include <iostream>

class InputInstance
{
public:
    InputInstance();
    
    friend std::ostream& operator<<(std::ostream& out, const InputInstance& instance);
    friend std::istream& operator>>(std::istream& in, InputInstance& instance);
    
    const Int3Array& C() const
    {
        return _C;
    }
    
    int k() const
    {
        return _k;
    }
    
    int numChr() const
    {
        return _num_chr;
    }
    
    const IntArray& n() const
    {
        return _n;
    }
    
    int e() const
    {
        int max_e = 0;
        for (int chr = 0; chr < _num_chr; ++chr)
        {
            for (int i = 0; i < _k; ++i)
            {
                for (int s = 0; s < _n[chr]; ++s)
                {
                    if (_C[chr][i][s] > max_e)
                    {
                        max_e = _C[chr][i][s];
                    }
                }
            }
        }
        return max_e;
    }
    
private:
    Int3Array _C;
    int _k;
    int _num_chr;
    IntArray _n;
};

std::ostream& operator<<(std::ostream& out, const InputInstance& instance);
std::istream& operator>>(std::istream& in, InputInstance& instance);

#endif // _INPUTINSTANCE_H