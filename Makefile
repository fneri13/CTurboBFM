# =========================
# Directories
# =========================
SRC_DIR := src
INC_DIR := include
OBJ_DIR := build
BIN_DIR := bin

# =========================
# External library paths
# =========================
# CoolProp
COOLPROP_INCLUDE_PATH := /Users/fneri/Documents/CoolProp/include
COOLPROP_LIB_PATH := /Users/fneri/Documents/CoolProp/build  # for static library

# fmt
FMT_INCLUDE_PATH := /opt/homebrew/include
FMT_LIB_PATH     := /opt/homebrew/lib

# gperftools (for profiling)
PROFILE_LIB_PATH := /opt/homebrew/opt/gperftools/lib
PROFILE_INCLUDE_PATH := /opt/homebrew/opt/gperftools/include

# =========================
# Compiler and flags
# =========================
CXX := g++-14
CXXFLAGS := -std=c++20 -Wall \
            -I$(INC_DIR) -I$(SRC_DIR) \
            -I$(COOLPROP_INCLUDE_PATH) \
            -I$(FMT_INCLUDE_PATH)
LDFLAGS := -L$(FMT_LIB_PATH) -lfmt

# Add CoolProp static library linking if you are using the compiled .a version
USE_COOLPROP_STATIC := 1
ifeq ($(USE_COOLPROP_STATIC),1)
    LDFLAGS += -L$(COOLPROP_LIB_PATH) -lCoolProp
endif

# Build type flags
DEBUG_FLAGS := -g -O0
RELEASE_FLAGS := -O3 -DNDEBUG -march=native
PROFILE_FLAGS := -g -O3 -DNDEBUG -march=native -I$(PROFILE_INCLUDE_PATH)
PROFILE_LDFLAGS := -L$(PROFILE_LIB_PATH) -lprofiler

# Default build type
BUILD_TYPE := debug

# =========================
# Sources and objects
# =========================
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
TARGET := $(BIN_DIR)/CTurboBFM

# =========================
# Build rules
# =========================
all: $(BUILD_TYPE)

debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(TARGET)

release: CXXFLAGS += $(RELEASE_FLAGS)
release: $(TARGET)

profile: CXXFLAGS += $(PROFILE_FLAGS)
profile: LDFLAGS += $(PROFILE_LDFLAGS)
profile: $(TARGET)

# =========================
# Ensure directories exist
# =========================
$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link object files into final binary
$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(OBJECTS) $(CXXFLAGS) $(LDFLAGS) -o $@

# =========================
# Clean
# =========================
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# =========================
# Phony targets
# =========================
.PHONY: all debug release profile clean
