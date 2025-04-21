rm -f build
mkdir build

clang++ -std=c++20 -g -Wall -c main.cpp -o build/main.o
clang++ -std=c++20 -g -Wall -c ../../src/CMesh.cpp -o build/CMesh.o
clang++ -std=c++20 -g -Wall -c ../../src/Config.cpp -o build/Config.o
clang++ -std=c++20 -g -Wall build/main.o build/CMesh.o build/Config.o -o build/main.out


