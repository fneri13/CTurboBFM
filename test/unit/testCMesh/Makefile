# Makefile

CXX = g++
CXXFLAGS = -std=c++20 -I ~/libraries/include/ -g -Wall
LDFLAGS = -L ~/libraries/lib -lgtest

SRC_DIR = ../../../src
BUILD_DIR = build
INCLUDES = ~/libraries/include/

OBJS = $(BUILD_DIR)/main.o $(BUILD_DIR)/CMesh.o $(BUILD_DIR)/Config.o $(BUILD_DIR)/commonFunctions.o
TARGET = $(BUILD_DIR)/main.out

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/main.o: main.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/CMesh.o: $(SRC_DIR)/CMesh.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/Config.o: $(SRC_DIR)/Config.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/commonFunctions.o: $(SRC_DIR)/commonFunctions.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR)
