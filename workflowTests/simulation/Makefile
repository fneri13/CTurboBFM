CXX = g++
SRC_DIR = ../../src
INCLUDES = ~/libraries/include/
VTK_INCLUDE = /opt/homebrew/include/vtk-9.4  # Adjust this path if necessary
VTK_LIB_DIR = /opt/homebrew/lib  # Adjust this path if necessary

# VTK Libraries to link (modify this if you need more libraries)
VTK_LIBS = -lvtkCommonCore-9.4 -lvtkIOXML-9.4 -lvtkCommonDataModel-9.4 -lvtkCommonExecutionModel-9.4 -lvtkIOCore-9.4 -lvtkCommonSystem-9.4 -lvtkCommonMath-9.4

# Source files
SOURCES = \
	main.cpp \
	$(SRC_DIR)/CMesh.cpp \
	$(SRC_DIR)/Config.cpp \
	$(SRC_DIR)/commonFunctions.cpp \
	$(SRC_DIR)/CSolverBase.cpp \
	$(SRC_DIR)/CEulerSolver.cpp \
	$(SRC_DIR)/CFluid.cpp \
	$(SRC_DIR)/CJSTScheme.cpp \
	$(SRC_DIR)/CBoundaryConditionEulerWall.cpp \
	$(SRC_DIR)/CBoundaryConditionInlet.cpp \
	$(SRC_DIR)/CBoundaryConditionOutlet.cpp \

# Common compiler flags
INCLUDE_FLAGS = -I $(INCLUDES) -I $(VTK_INCLUDE)

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

# Link final binaries with VTK libraries
$(DEBUG_TARGET): $(DEBUG_OBJS)
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_OBJS) $(VTK_LIB_DIR)/libvtkCommonCore-9.4.a $(VTK_LIB_DIR)/libvtkIOXML-9.4.a $(VTK_LIB_DIR)/libvtkCommonDataModel-9.4.a $(VTK_LIB_DIR)/libvtkCommonExecutionModel-9.4.a $(VTK_LIB_DIR)/libvtkIOCore-9.4.a $(VTK_LIB_DIR)/libvtkCommonSystem-9.4.a $(VTK_LIB_DIR)/libvtkCommonMath-9.4.a -o $@

$(RELEASE_TARGET): $(RELEASE_OBJS)
	$(CXX) $(RELEASE_FLAGS) $(RELEASE_OBJS) $(VTK_LIB_DIR)/libvtkCommonCore-9.4.a $(VTK_LIB_DIR)/libvtkIOXML-9.4.a $(VTK_LIB_DIR)/libvtkCommonDataModel-9.4.a $(VTK_LIB_DIR)/libvtkCommonExecutionModel-9.4.a $(VTK_LIB_DIR)/libvtkIOCore-9.4.a $(VTK_LIB_DIR)/libvtkCommonSystem-9.4.a $(VTK_LIB_DIR)/libvtkCommonMath-9.4.a -o $@

clean:
	rm -rf build
