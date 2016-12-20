# ILP for the Copy-Number Tree (CNT) Problem
============================================

## Dependencies

* [LEMON Graph library](http://lemon.cs.elte.hu) (>= 1.3)
* [ILOG CPLEX](http://www.ibm.com/developerworks/downloads/ws/ilogcplex/) (>= 12.0)

## Compilation instructions

To compile `CNT`, execute the following commands from the root of the repository:

    mkdir build
    cd build
    cmake ..
    make
    
Note: On Mac OS >= 10.9, you have to ensure that LEMON is linked against the old C++ standard library by editing LEMON's `CMakeLists.txt` file as follows.

	if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
	  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
	endif()
    
In case CMake fails to detect either CPLEX or LEMON, run the following command with adjusted paths:

	cmake \
	-DLIBLEMON_ROOT=~/lemon \
	-DCPLEX_INC_DIR=~/ILOG/cplex/include/ \
	-DCPLEX_LIB_DIR=~/ILOG/cplex/lib/x86-64_osx/static_pic \
	-DCONCERT_LIB_DIR=~/ILOG/concert/lib/x86-64_osx/static_pic \
	-DCONCERT_INC_DIR=~/ILOG/concert/include/ ..
	
The compilation results in the following files in `build` directory:

* `cnt`
* `visualize`

## Usage instructions

### Input format

Example:

    #PARAMS
    1 #number of chromosomes
    4 #number of leaves
    10 #number of segments for each chromosome
    #PROFILES
    3 : 2 2 3 3 3 3 3 3 2 2
    4 : 2 2 2 2 3 3 2 2 3 2
    5 : 1 1 1 1 2 2 1 1 1 1
    6 : 2 2 2 3 4 4 3 3 3 3

The first line must be `#PARAMS`. The three subsequent lines indicate the number of chromosomes, the number *k* of leaves and the number of segments for each chromosome (separated by a space). The fourth line must be `#PROFILES`. Then for each leaf the copy number profile is given (separated by spaces). It is recommended to label the leaves by the numbers k-1, ..., 2k-2.

### Example executions

To run `cnt`:

    ./cnt -s 2 ../data/20160430/simC_c1_k4_n10_u1_s1_del02.input > result.txt

To visualize the resulting copy-number tree:

    ./visualize result.txt > T.dot
    dot -Tpdf T.dot -o T.pdf

