include make.inc
# Group definitions
default: test.exe

# Executable compilations
test.exe: test.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

# Object file creations

test.o: test.c
	$(CC) -c $(CFLAGS) $^ -o $@
