include make.inc
# Group definitions
default: test.exe test_optBlas.exe

clean:
	rm *o driver *exe

# Executable compilations
test_optBlas.exe: test.o test.c
	$(FC) $(CFLAGS) test.o $(OPTBLAS) -o $@

test.exe: test.o test.c
	$(FC) $(CFLAGS) test.o $(LAPACKELIB) $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

# Object file creations
test.o: test.c 
	$(CC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^ -o $@
