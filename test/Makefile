
# Directories
SRC_DIR = ../src
INC_DIR = ../include
LIB_DIR = ../lib

# Set compiler variables
CXX = mpic++
CXXFLAGS = -std=c++20 -Wall -fopenmp
CPPFLAGS = -I$(INC_DIR) -I$(INC_DIR)/muparser

# Set linker variables
LDFLAGS ?= -L$(LIB_DIR) -fopenmp
LIBS ?= -lmuparser -Wl,-rpath=$(LIB_DIR)

# Get all source files in the src directory
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:.cpp=.o)
HEADERS = $(wildcard $(INC_DIR)/*.hpp)

# Test files
TEST_SRCS = $(wildcard *.cpp)
TEST_OBJS = $(TEST_SRCS:.cpp=.o)
EXEC = $(TEST_SRCS:.cpp=)

.phony= all clean distclean doc no_omp

.DEFAULT_GOAL = all

all: $(OBJS)

test: $(EXEC)

# Build the executables
%: %.o $(OBJS)
$(EXEC) : $(TEST_OBJS) $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) $(LIBS) -o $@

# Compile objects files from src
$(SRC_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -c -o $@
	
# Compile object files from test
%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Recompile without _OPENMP
no_omp:
	$(MAKE) clean
	$(MAKE) test CXXFLAGS="$(CXXFLAGS:-fopenmp=)" LDFLAGS="$(LDFLAGS:-fopenmp=)"


clean:
	$(RM) -f $(OBJS) $(TEST_OBJS)

distclean: clean
	$(RM) $(EXEC)
	$(RM) *.out *.vtk performance_results.txt *~ -r ../doc
	$(RM) -r ../doc

doc:
	doxygen $(DOXYFILE)
	$(MAKE) -C ../doc/latex
