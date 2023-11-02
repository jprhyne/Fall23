include make.inc

default: test.exe

createMat.o: createMat.f
	$(FC) -c $(FFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^

my_dorgqr.o: my_dorgqr_v1.f
	$(FC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^ -o $@

test.o: test.c 
	$(CC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^

test.exe: my_dorgqr.o test.o createMat.o test.c
	$(FC) $(CFLAGS) my_dorgqr.o test.o createMat.o $(OPTBLAS) -o $@

clean:
	rm *o test.exe
