CXX = g++
CXXFLAGS = -std=c++17 -Wall -I ../../../include -I ~/libraries/include
LDFLAGS = -L ~/libraries/lib
LDLIBS = -lgtest -lgtest_main -pthread

SRC = ../../../src/input.cpp ../../../src/mesh.cpp ../../../src/config.cpp ../../../src/math_utils.cpp main.cpp
OBJ = $(SRC:.cpp=.o)

TARGET = main

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)