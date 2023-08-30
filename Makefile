include make.inc

default: test.exe

my_dorgqr.o: my_dorgqr.f
	$(FC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^

test.o: test.c
	$(CC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^

test.exe: my_dorgqr.o test.o test.c
	$(FC) $(CFLAGS) my_dorgqr.o test.o $(LAPACKELIB) $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@
