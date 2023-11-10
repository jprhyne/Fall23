include make.inc

default: test.exe

null.o: null.f
	$(FC) -c $(FFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^

createMat.o: createMat.f
	$(FC) -c $(FFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^

computeQ2.o: computeQ2.f
	$(FC) -c $(FFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^

my_dorgqr.o: my_dorgqr_v2.f
	$(FC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^ -o $@

test.o: test.c 
	$(CC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^

test.exe: my_dorgqr.o test.o createMat.o computeQ2.o test.c null.o
	$(FC) $(CFLAGS) my_dorgqr.o test.o createMat.o computeQ2.o null.o $(OPTBLAS) -o $@

clean:
	rm *o test.exe
