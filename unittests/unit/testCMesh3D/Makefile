# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -g -Wall -I../../../include/ -I$(HOME)/libraries/include/
LDFLAGS = -L$(HOME)/libraries/lib -lgtest

# Directories
SRC_DIR = ../../../src
BUILD_DIR = build

# Object files
OBJS = $(BUILD_DIR)/main.o \
       $(BUILD_DIR)/CMesh.o \
       $(BUILD_DIR)/Config.o \
       $(BUILD_DIR)/commonFunctions.o \
       $(BUILD_DIR)/CFluid.o

# Output binary
TARGET = $(BUILD_DIR)/main.out

# Default target
all: $(TARGET)

# Ensure build directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compilation rules
$(BUILD_DIR)/main.o: main.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/CMesh.o: $(SRC_DIR)/CMesh.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/Config.o: $(SRC_DIR)/Config.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/commonFunctions.o: $(SRC_DIR)/commonFunctions.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/CFluid.o: $(SRC_DIR)/CFluid.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Linking
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Clean rule
clean:
	rm -rf $(BUILD_DIR)
