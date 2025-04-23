# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -g -Wall

# Source files
SRCS = main.cpp ../../src/CMesh.cpp ../../src/Config.cpp ../../src/commonFunctions.cpp ../../src/CSolverBase.cpp ../../src/CEulerSolver.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Output folder
BUILD_DIR = build

# Final executable
EXEC = $(BUILD_DIR)/main.out

# Default target
all: $(BUILD_DIR) $(EXEC)

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Link object files to create the final executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC)

# Compile .cpp files to .o files
$(BUILD_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up the build folder
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
