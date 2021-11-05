# esc412_advanced_high_performance_computing

Code written by me can be found in src/.


## Build the project

First build blitz.
```bash
cmake -DCMAKE_INSTALL_PREFIX=./ -B ./blitz-build ./blitz-1.0.2
make -C ./blitz-build
make -C ./blitz-build install
```

Create the necessary folder
```bash
mkdir input
mkdir output
mkdir run
```
Put the particle file in input.

Then either run the project with clion and the provide CMakeLists.txt in the root folder or
compile and run from the root folder in the command line with `g++` or other c++ compiler. 
```bash
 g++ -Wall -O2 -I include -std=c++17 -o run/main src/main.cpp src/tipsy.cpp src/my_io.cpp && ./run/main 
```

You might have to adjust the src files.