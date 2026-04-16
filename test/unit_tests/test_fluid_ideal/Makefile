# Compiler
CXX = g++
CXXFLAGS = -std=c++17 -Wall -I ../../../include/ -I ~/libraries/include

# Libraries
LDFLAGS = -L ~/libraries/lib
LDLIBS = -lgtest -lgtest_main -pthread

# Sources and objects
SRC = main.cpp ../../../src/fluid_ideal.cpp
OBJ = $(SRC:.cpp=.o)

# Target
TARGET = main

# Default rule
all: $(TARGET)

# Link step
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LDLIBS) -o $(TARGET)

# Compile step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJ) $(TARGET)