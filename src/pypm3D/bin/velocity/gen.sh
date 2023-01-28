# Generate lib
gcc -fPIC -Ofast -c velocity.c
gcc -shared -Ofast -o velocity.so velocity.o

# Remove intermediate file
rm velocity.o