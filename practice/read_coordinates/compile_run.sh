g++ -std=c++20 -c main.cpp -o main.o
g++ -std=c++20 -c ../../src/CMesh.cpp -o Cmesh.o

g++ -std=c++20 main.o Cmesh.o -o main.a

./main.a
