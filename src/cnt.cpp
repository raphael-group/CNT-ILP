#include <fstream>
#include <ios>
#include <stdlib.h>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include <lemon/arg_parser.h>

#include "basic_types.h"
#include "architect.h"
#include "rungarchitect.h"
#include "copynumbertree.h"

//void readK(const string &fileK,
//           const unsigned int num_chromosomes,
//           const IntArray &num_segments,
//           IntMatrix &K);

//options_t parse_arguments(const unsigned int argc, char** argv)
//{
//    options_t options;
//    
//    if (argc == 3)
//    {
//        options.options_initialized = true;
//        options.input_filename = argv[1];
//        int num = 0;
//        std::istringstream iss(argv[2]);
//        if((iss >> num).fail())
//        {
//            options.fileK = argv[2];
//        }
//        else {
//            options.K = num;
//        }
//    }
//    else {
//        std::cerr << "Usage: " << argv[0] << " ARG1 ARG2 where " << std::endl;
//        std::cerr << "ARG1: name of the input C file" << std::endl;
//        std::cerr << "ARG2: it can be equal to:" << std::endl;
//        std::cerr << "  1. An integer specifying the same maximum copy-number for any segment" << std::endl;
//        std::cerr << "  2. A string specifying the file name containing the copy-number for any segment" << std::endl;
//    }
//    
//    return options;
//}

int main(int argc, char** argv)
{
    int timeLimit = -1;
    int memoryLimit = -1;
    int maxCopyNumber = -1;
    bool rootNotFixed = false;
    int solver = 1;
    std::string exportFilename;

    lemon::ArgParser ap(argc, argv);
    ap.refOption("t", "Time limit in seconds (default: -1, disabled)", timeLimit)
      .refOption("m", "Memory limit in MB (default: -1, disabled)", memoryLimit)
      .refOption("e", "Maximum copy number (default: -1, inferred from leaves)", maxCopyNumber, false)
      .refOption("x", "Export LP filename", exportFilename)
      .refOption("r", "Do not fix root to all 2s", rootNotFixed)
      .refOption("s", "Choose the solver: (1) fake and real solver (2) rung solver", solver)
      .other("input", "Input file");
    ap.parse();
    
    if (ap.files().size() == 0)
    {
        std::cerr << "ERROR: missing input file" << std::endl;
        return 1;
    }
    
    InputInstance inputInstance;
    std::ifstream inFile(ap.files()[0].c_str());
    
    if (!inFile.good())
    {
        std::cerr << "ERROR: could not open '" << ap.files()[0] << "' for reading" << std::endl;
        return EXIT_FAILURE;
    }
    
    inFile >> inputInstance;
    
    if (maxCopyNumber < 0)
    {
        maxCopyNumber = inputInstance.e();
    }

    IntMatrix e(inputInstance.numChr());
    for(unsigned int chr = 0; chr < inputInstance.numChr(); ++chr)
    {
        e[chr] = IntArray(inputInstance.n()[chr], maxCopyNumber);
    }

    BaseArchitect* pArchitect = NULL;
    switch (solver)
    {
    case 1:
        pArchitect = new Architect(inputInstance, e, rootNotFixed);
        break;
    case 2:
        pArchitect = new RungArchitect(inputInstance, e, rootNotFixed);
        break;
    default:
        std::cerr << "Error: invalid solver specified" << std::endl;
        break;
    }

    pArchitect->init();
    
    if (!exportFilename.empty())
    {
        pArchitect->exportModel(exportFilename);
    }
    
    if (pArchitect->solve(timeLimit, memoryLimit))
    {
        std::cout << pArchitect->getTree();
    }

    delete pArchitect;
    return 0;
}


//void readK(const string &fileK,
//           const unsigned int num_chromosomes,
//           const IntArray &num_segments,
//           IntMatrix &e)
//{
//    std::ifstream input;
//    e = IntMatrix(num_chromosomes);
//    for(unsigned int chr = 0; chr < num_chromosomes; ++chr)
//    {
//        e[chr] = IntArray(num_segments[chr]);
//    }
//
//    try
//    {
//        input.open(fileK.c_str(), ios::in);
//    }
//    catch(const std::exception & ex)
//    {
//        std::cerr << "ERROR: failing opening the input file: " << fileK << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//    try
//    {
//        std::string chr;
//
//        unsigned int count_chr = 0;
//        while(!input.eof())
//        {
//            getline(input, chr, '\n');
//
//            stringstream schr(chr);
//            string value;
//            unsigned int count_value = 0;
//
//            while(!schr.eof()) {
//                std::getline(schr, value, ' ');
//                if(!value.empty())
//                {
//                    if(count_chr >= num_chromosomes || count_value >= num_segments[count_chr])
//                    {
//                        std::cerr << "ERROR: : the K format of " << fileK
//                                  << "is wrong or the number of chromosomes and corresponding segments is inconsistent with the input" << std::endl;
//                        exit(EXIT_FAILURE);
//                    }
//                    else
//                    {
//                        if(atof(value.c_str()) < 2)
//                        {
//                            throw std::runtime_error("ERROR: K must be at least 2");
//                            exit(EXIT_FAILURE);
//                        }
//                        else
//                        {
//                            e[count_chr][count_value] = atof(value.c_str());
//                        }
//                    }
//                    ++count_value;
//                }
//            }
//            ++count_chr;
//        }
//    }
//    catch(const std::exception & ex)
//    {
//        std::cerr << "ERROR: failing opening the input file: " << fileK << "\": " << ex.what() << std::endl;
//        exit(EXIT_FAILURE);
//    }
//}
