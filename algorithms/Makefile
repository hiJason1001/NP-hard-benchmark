CXX = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -Wreorder -Wunused -Wshadow -Werror=return-type -std=c++20

# SRC = $(wildcard *.cpp)
SRC = NN.cpp

# Derive output file names based on the source file names
OUT = $(SRC:.cpp=.out)

all: $(OUT)

# create an output file from a source file
%.out: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(OUT)

rebuild: clean all
