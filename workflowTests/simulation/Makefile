CXX = g++
SRC_DIR = ../../src
INCLUDES = ~/libraries/include/

# Source files
SOURCES = \
	main.cpp \
	$(SRC_DIR)/CMesh.cpp \
	$(SRC_DIR)/Config.cpp \
	$(SRC_DIR)/commonFunctions.cpp \
	$(SRC_DIR)/CSolverBase.cpp \
	$(SRC_DIR)/CEulerSolver.cpp \
	$(SRC_DIR)/CFluid.cpp \
	$(SRC_DIR)/CJSTScheme.cpp

# Common compiler flags
INCLUDE_FLAGS = -I $(INCLUDES)

# Build types
DEBUG_DIR = build/debug
RELEASE_DIR = build/release

DEBUG_FLAGS = -std=c++20 $(INCLUDE_FLAGS) -g -O0 -Wall
RELEASE_FLAGS = -std=c++20 $(INCLUDE_FLAGS) -O3 -Wall

DEBUG_OBJS = $(patsubst %.cpp, $(DEBUG_DIR)/%.o, $(notdir $(SOURCES)))
RELEASE_OBJS = $(patsubst %.cpp, $(RELEASE_DIR)/%.o, $(notdir $(SOURCES)))

DEBUG_TARGET = $(DEBUG_DIR)/main.out
RELEASE_TARGET = $(RELEASE_DIR)/main.out

.PHONY: all debug release clean

all: release

debug: $(DEBUG_TARGET)
release: $(RELEASE_TARGET)

# Make sure the folders exist
$(DEBUG_DIR):
	mkdir -p $(DEBUG_DIR)

$(RELEASE_DIR):
	mkdir -p $(RELEASE_DIR)

# Compile each .cpp to .o in the debug folder
$(DEBUG_DIR)/%.o: $(SRC_DIR)/%.cpp | $(DEBUG_DIR)
	$(CXX) $(DEBUG_FLAGS) -c $< -o $@

$(DEBUG_DIR)/main.o: main.cpp | $(DEBUG_DIR)
	$(CXX) $(DEBUG_FLAGS) -c $< -o $@

# Compile each .cpp to .o in the release folder
$(RELEASE_DIR)/%.o: $(SRC_DIR)/%.cpp | $(RELEASE_DIR)
	$(CXX) $(RELEASE_FLAGS) -c $< -o $@

$(RELEASE_DIR)/main.o: main.cpp | $(RELEASE_DIR)
	$(CXX) $(RELEASE_FLAGS) -c $< -o $@

# Link final binaries
$(DEBUG_TARGET): $(DEBUG_OBJS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_OBJS) -o $@

$(RELEASE_TARGET): $(RELEASE_OBJS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_OBJS) -o $@

clean:
	rm -rf build
