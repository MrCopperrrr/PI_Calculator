## Build
g++ -O3 -march=native main.cpp -o main -lmpfr -lgmpxx -lgmp
g++ -O3 -march=native -flto  main.cpp -o main -lmpfr -lgmpxx -lgmp
g++ -Ofast -march=native -mtune=native -funroll-loops -flto -m64 -fopenmp -fomit-frame-pointer -DNDEBUG -pipe -o main main.cpp -lmpfr -lgmpxx -lgmp -pthread


## Run
./main 16 100000