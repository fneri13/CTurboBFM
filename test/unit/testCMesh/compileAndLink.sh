# Clean and create the build folder
rm -r build
mkdir build

# Compile C++ files
g++ -std=c++20 -I ~/libraries/include/ -g -Wall -c main.cpp -o build/main.o
g++ -std=c++20 -I ~/libraries/include/ -g -Wall -c ../../../src/CMesh.cpp -o build/CMesh.o
g++ -std=c++20 -I ~/libraries/include/ -g -Wall -c ../../../src/Config.cpp -o build/Config.o

# Link object files to create the final executable
g++ -std=c++20 -g -Wall build/main.o build/CMesh.o build/Config.o -o build/main.out -L ~/libraries/lib -lgtest


