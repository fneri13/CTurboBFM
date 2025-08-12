# Directories
SRC_DIR := src
INC_DIR := include
OBJ_DIR := build
BIN_DIR := bin

PROFILE_LIB_PATH := /opt/homebrew/opt/gperftools/lib
PROFILE_INCLUDE_PATH := /opt/homebrew/opt/gperftools/include

# Compiler and flags
CXX := g++-14
CXXFLAGS := -std=c++20 -Wall -I$(INC_DIR) -I$(SRC_DIR)

# Debug and Release specific flags
DEBUG_FLAGS := -g -O0
RELEASE_FLAGS := -O3 -DNDEBUG -march=native
PROFILE_FLAGS := -g -O3 -DNDEBUG -march=native -I$(PROFILE_INCLUDE_PATH)

# Default build type
BUILD_TYPE := debug

# Sources and derived object file names
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))

# Output binary
TARGET := $(BIN_DIR)/CTurboBFM

# Default rule
all: $(TARGET)

debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(TARGET)

release: CXXFLAGS += $(RELEASE_FLAGS)
release: $(TARGET)

profile: CXXFLAGS += $(PROFILE_FLAGS)
profile: LDFLAGS += -L$(PROFILE_LIB_PATH) -lprofiler
profile: $(TARGET)

# Link object files into final binary
$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -o $(TARGET)

# Compile source files to object files in build/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean debug release install
