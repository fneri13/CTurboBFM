g++ -std=c++20 -Wall -g -c main.cpp -o main.o
g++ -std=c++20 -Wall -g *.o -o main.a

rm *.o
./main.a