# Finite Size Scaling v.0.1.0

## create a C file:


gcc -fPIC -shared -o src/simulation.so src/all.c -I/home/emilio/Projects/finite_size_scaling/src/aux/ -I/user/lib64/include/ -lgsl -lgslcblas -lm -L/user/lib64/lib -O2