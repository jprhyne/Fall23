include make.inc
# Group definitions
default: test.exe test_optBlas.exe

clean:
	rm *o *exe

# Executable compilations
test_optBlas.exe: test_opt.o test.c
	$(FC) -D USE_AOCL $(CFLAGS) test_opt.o $(OPTBLAS) -o $@

test.exe: test.o test.c
	$(FC) $(CFLAGS) test.o $(LAPACKELIB) $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

# Object file creations
test.o: test.c 
	$(CC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^ -o $@

test_opt.o: test.c 
	$(CC) -D USE_AOCL -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^ -o $@
