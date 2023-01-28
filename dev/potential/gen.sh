# Generate lib
gcc -fPIC -Ofast -c potential.c
gcc -shared -Ofast -o potential.so potential.o

# Remove intermediate file
rm potential.o