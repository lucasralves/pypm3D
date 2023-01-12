# Generate lib
gcc -fPIC -Ofast -c lib.c
gcc -shared -Ofast -o lib.so lib.o

# Remove intermediate file
rm lib.o