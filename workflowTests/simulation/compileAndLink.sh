# Clean and create the build folder
rm -r build
mkdir build

# Compile C++ files
g++ -std=c++20 -g -Wall -c main.cpp -o build/main.o
g++ -std=c++20 -g -Wall -c ../../src/CMesh.cpp -o build/CMesh.o
g++ -std=c++20 -g -Wall -c ../../src/Config.cpp -o build/Config.o
g++ -std=c++20 -g -Wall -c ../../src/commonFunctions.cpp -o build/commonFunctions.o

# Link object files to create the final executable
g++ -std=c++20 -g -Wall build/*.o -o build/main.out


