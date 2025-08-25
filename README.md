## Build dependencies (Ubuntu/Debian)
```sh
sudo apt update
sudo apt install build-essential libgmp-dev libmpfr-dev
```

## Build
g++ -O3 -march=native main.cpp -o main -lmpfr -lgmpxx -lgmp


## Run
./main 16 100000

## Note
This program can only calculate up to 100 million digits of pi.