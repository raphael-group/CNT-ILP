#include "inputinstance.h"
#include <sstream>
#include <string>
#include <stdexcept>

InputInstance::InputInstance()
    : _C()
    , _k(-1)
    , _num_chr(-1)
{
}

std::ostream& operator<<(std::ostream& out, const InputInstance& instance)
{
    out << "#PARAMS" << std::endl;
    out << instance._num_chr << " #number of chromosomes" << std::endl;
    out << instance._k << " #number of leaves" << std::endl;
    for (int n : instance._n)
    {
        out << n << " ";
    }
    out << "#number of segments per chromosome" << std::endl;
    out << "#PROFILES" << std::endl;

    for (int i = 0; i < instance._k; ++i)
    {
        out << i << " :";
        for (int chr = 0; chr < instance._num_chr; ++chr)
        {
            if (chr != 0)
            {
                out << " |";
            }
            for (int s = 0; s < instance._n[chr]; ++s)
            {
                out << " " << instance._C[chr][i][s];
            }
        }
        out << std::endl;
    }
    return out;
}

std::istream& operator>>(std::istream& in, InputInstance& instance)
{
    std::string line;
    std::string value;
    
    instance._C.clear();
    instance._n.clear();
    
    std::getline(in, line, '\n'); //Skip the first line "#PARAMS"
    
    /** Read the number of samples **/
    std::getline(in, line, '\n');
    std::stringstream chromosomeline(line);
    chromosomeline >> instance._num_chr;
    
    /** Read the number of chromosomes **/
    std::getline(in, line, '\n');
    std::stringstream leafline(line);
    leafline >> instance._k;
    
    /** Read the number of segments for each chromosome  **/
    std::getline(in, line, '\n');
    std::stringstream segline(line);
    std::getline(segline, line, '#');
    std::stringstream split_seg(line);
    while (!split_seg.eof())
    {
        std::getline(split_seg, value, ' ');
        if(!value.empty())
        {
            instance._n.push_back((unsigned int)atoi(value.c_str()));
        }
    }
    
    if(instance._n.size() != instance._num_chr)
    {
        throw std::runtime_error(std::string("ERROR: inconsistent numbers of segments in #PARAMS"));
    }
    
    /** Allocate the cube containing the result **/
    instance._C = Int3Array(instance._num_chr, IntMatrix(instance._k));
    
    for(unsigned int chr = 0; chr < instance._num_chr; ++chr)
    {
        for(unsigned int leaf = 0; leaf < instance._k; ++leaf)
        {
            instance._C[chr][leaf] = IntArray(instance._n[chr]);
        }
    }
    
    getline(in, line, '\n'); //Skip the line "#PROFILES"
    
    /** Begin the reading of the profiles **/
    unsigned int counter_chromosomes = 0;
    unsigned int counter_leaves = 0;
    unsigned int counter_seg = 0;
    
    while (counter_leaves < instance._k && in.good())
    {
        getline(in, line, '\n');
        if(!line.empty() && !(line == "#EDGES"))
        {
            std::stringstream csline(line);
            getline(csline, line, ':'); //Skip the name Li of the leaf
            getline(csline, line, '\n');
            std::stringstream sline(line);
            
            counter_chromosomes = 0;
            while(!sline.eof())
            {
                std::string leaf;
                std::getline(sline, leaf, '|');
                std::stringstream sleaf(leaf);
                
                counter_seg = 0;
                while(!sleaf.eof())
                {
                    getline(sleaf, value, ' ');
                    if(!value.empty())
                    {
                        if(counter_leaves >= instance._k
                           || counter_chromosomes >= instance._num_chr
                           || counter_seg >= instance._n[counter_chromosomes])
                        {
                            throw std::runtime_error("ERROR: the input C format is wrong or the number of chromosomes/leaves/corresponding segments is inconsistent with the specified PARAMS");
                        }
                        else
                        {
                            instance._C[counter_chromosomes][counter_leaves][counter_seg] = atoi(value.c_str());
                            
                        }
                        ++counter_seg;
                    }
                }
                
                if(counter_seg != instance._n[counter_chromosomes])
                {
                    std::stringstream tmp;
                    tmp << "ERROR: inconsistent number of segments for chromosome " << counter_chromosomes;
                    
                    throw std::runtime_error(tmp.str());
                }
                
                ++counter_chromosomes;
            }
            
            if(counter_chromosomes != instance._num_chr)
            {
                throw std::runtime_error("ERROR: inconsistent number of chromosomes");
            }
            
            ++counter_leaves;
        }
    }
    
    if(counter_leaves != instance._k)
    {
        throw std::runtime_error("ERROR: inconsistent number of leaves");
    }

    return in;
}
