rm *.a 

g++ -std=c++20 -c -Wall -g main.cpp -o main.o
g++ -std=c++20 main.o -o main.a

rm *.o
./main.a
