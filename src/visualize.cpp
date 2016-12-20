#include "basic_types.h"
#include "copynumbertree.h"
#include <string>
#include <fstream>

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <FILE> where" << std::endl
                  << "  <FILE> is copy-number tree filename or - for STDIN" << std::endl;
        
        return 1;
    }
    
    std::string filename(argv[1]);

    CopyNumberTree T;
    if (filename == "-")
    {
        std::cin >> T;
    }
    else
    {
        std::ifstream inFile(filename.c_str());
        if (!inFile.good())
        {
            std::cerr << "ERROR: could not open '" << filename << "' for reading" << std::endl;
            return 1;
        }
        inFile >> T;
        inFile.close();
    }
    
    T.writeDOT(std::cout);
    
    return 0;
}
