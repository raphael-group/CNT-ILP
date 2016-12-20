#ifndef _BASIC_H_
#define _BASIC_H_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <set>
#include <algorithm>

typedef std::vector<unsigned int> IntArray;
typedef std::vector<IntArray> IntMatrix;
typedef std::vector<IntMatrix> Int3Array;
typedef std::vector<Int3Array> Int4Array;

typedef std::vector<double> DoubleArray;
typedef std::vector<DoubleArray> DoubleMatrix;
typedef std::vector<DoubleMatrix> Double3Array;
typedef std::vector<Double3Array> Double4Array;


typedef std::set<int> IntSet;
typedef std::pair<IntSet, IntSet> IntSetPair;


//struct options_t
//{
//    bool options_initialized;
//    std::string input_filename;
//    std::string fileK;
//    unsigned int K;
//
//    options_t()
//    : options_initialized(false)
//    , input_filename("")
//    , fileK("")
//    , K(0)
//    {}
//};

//inline
//std::string itos(const int i)
//{
//  std::ostringstream result;
//  result << i;
//  return result.str();
//}

#endif // _BASIC_H_
