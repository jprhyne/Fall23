include make.inc
# Group definitions
default: test.exe driver test_optBlas.exe timeSVD.exe

testDorgqr: test_mklDorgqr.exe test_refDorgqr.exe

nanTest: test_v1.exe test_v1_opt.exe

clean:
	rm *o driver *exe

# Executable compilations
timeSVD_opt.exe: timeSVD.o my_dorgqr_v1.o my_dgesvd.o timeSVD.c
	$(FC) $(CFLAGS) timeSVD.o my_dorgqr_v1.o my_dgesvd.o $(LAPACKELIB) $(OPTBLAS) -o $@

timeSVD.exe: timeSVD.o my_dorgqr.o my_dgesvd.o dgesvd.o timeSVD.c
	$(FC) $(CFLAGS) timeSVD.o my_dorgqr.o my_dgesvd.o dgesvd.o $(LAPACKLIB) $(BLASLIB) $(LAPACKELIB) -o $@

test_optBlas.exe: my_dorgqr.o test.o test.c
	$(FC) $(CFLAGS) my_dorgqr.o test.o $(OPTBLAS) -o $@

test_refDorgqr.exe: test.o dorgqr.o test.c
	$(FC) $(CFLAGS) test.o $(OPTBLAS) dorgqr.o -o $@

test_mklDorgqr.exe: test.o test.c
	$(FC) $(CFLAGS) test.o $(OPTBLAS) -o $@

test_v2.exe: test.o my_dorgqr_v2.o test.c
	$(FC) $(CFLAGS) test.o my_dorgqr_v2.o $(OPTBLAS) -o $@

test_v1_opt.exe: test.o my_dorgqr_v1.o test.c
	$(FC) $(CFLAGS) test.o my_dorgqr_v1.o $(OPTBLAS) -o $@

test_v1.exe: test.o my_dorgqr_v1.o test.c
	$(FC) $(CFLAGS) test.o my_dorgqr_v1.o $(LAPACKELIB) $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

test.exe: my_dorgqr.o test.o test.c
	$(FC) $(CFLAGS) my_dorgqr.o test.o $(LAPACKELIB) $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

driver: driver.c
	$(CC) $(CFLAGS) $^ -o $@

# Object file creations
timeSVD.o: timeSVD.c
	$(CC) -c $(CFLAGS) -DUSE_AMD=1 $^ -o $@

my_dgesvd.o: my_dgesvd.f
	$(FC) -c $(FFLAGS) $^ -o $@

dgesvd.o: dgesvd.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr.o: my_dorgqr.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v1.o: my_dorgqr_v1.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v2.o: my_dorgqr_v2.f
	$(FC) -c $(FFLAGS) $^ -o $@

test.o: test.c 
	$(CC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^ -o $@

dorgqr.o: dorgqr.f
	$(FC) -c $(CFLAGS) $^ -o $@
