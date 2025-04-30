# Directories
SRC_DIR := src
INC_DIR := include
OBJ_DIR := build
BIN_DIR := bin
INSTALL_DIR := /usr/local/bin  # Modify this as needed (e.g., use ~/bin or another directory)

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -g -I$(INC_DIR) -I$(SRC_DIR)

# Debug and Release specific flags
DEBUG_FLAGS := -g -O0
RELEASE_FLAGS := -O3 -DNDEBUG

# Default build type
BUILD_TYPE := debug

# Sources and derived object file names
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))

# Output binary
TARGET := $(BIN_DIR)/CTurboBFM

# Default rule
all: $(TARGET)

# Build Debug mode
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(TARGET)

# Build Release mode
release: CXXFLAGS += $(RELEASE_FLAGS)
release: $(TARGET)

# Link object files into final binary
$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJECTS) -o $(TARGET)

# Compile source files to object files in build/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Install the binary to the specified directory
install: $(TARGET)
	@mkdir -p $(INSTALL_DIR)
	cp $(TARGET) $(INSTALL_DIR)
	@echo "Installed $(TARGET) to $(INSTALL_DIR)"

.PHONY: all clean debug release install
