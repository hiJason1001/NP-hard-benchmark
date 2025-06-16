CXX = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -Wreorder -Wunused -Wshadow -Werror=return-type -std=c++20

# SRC = $(wildcard *.cpp)
SRC = HPA.cpp

OUT = $(SRC:.cpp=.out)

all: $(OUT)

%.out: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(OUT)

rebuild: clean all
