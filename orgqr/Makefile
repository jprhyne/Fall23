include make.inc
# Group definitions
default: test.exe driver test_optBlas.exe timeSVD.exe tests

tests: test_v0.exe test_v1_opt.exe test_v2.exe test_v3.exe test_v4.exe test_v5.exe test_v7.exe

testDorgqr: test_mklDorgqr.exe test_refDorgqr.exe

nanTest: test_v1.exe test_v1_opt.exe

clean:
	rm *o driver *exe

# Executable compilations
timeSVD_opt.exe: timeSVD.o my_dorgqr_v1.o my_dgesvd.o timeSVD.c
	$(FC) $(CFLAGS) timeSVD.o my_dorgqr_v1.o my_dgesvd.o $(LAPACKELIB) $(OPTBLAS) -o $@

timeSVD.exe: timeSVD.o my_dorgqr.o my_dgesvd.o dgesvd.o timeSVD.c
	$(FC) $(CFLAGS) timeSVD.o my_dorgqr.o my_dgesvd.o dgesvd.o $(LAPACKLIB) $(BLASLIB) $(LAPACKELIB) -o $@

test_optBlas.exe: my_dorgqr.o test_my_dorgqr.o
	$(FC) $(CFLAGS) $^ $(LAPACKELIB) $(OPTBLAS) -o $@

test_refDorgqr.exe: test.o
	$(FC) $(CFLAGS) test.o $(OPTBLAS) -o $@

test_mklDorgqr.exe: test.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

test_v0.exe: test_my_dorgqr.o my_dorgqr_v0.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

test_v1.exe: test_my_dorgqr.o my_dorgqr_v1.o
	$(FC) $(CFLAGS) $^ $(LAPACKELIB) $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

test_v1_opt.exe: test_my_dorgqr.o my_dorgqr_v1.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

test_v2.exe: test_my_dorgqr.o my_dorgqr_v2.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

test_v3.exe: test_my_dorgqr.o my_dorgqr_v3.o ../org2r/src/my_dorgkr.o ../org2r/src/lumm.o ../org2r/src/dtvt.o ../org2r/src/dtrmmoop.o ../org2r/src/my_dlarft.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

test_v4.exe: test_my_dorgqr.o my_dorgqr_v4.o ../org2r/src/my_dorgkr.o ../org2r/src/lumm.o ../org2r/src/dtvt.o ../org2r/src/dtrmmoop.o ../org2r/src/my_dlarft.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

test_v5.exe: test_my_dorgqr.o my_dorgqr_v5.o my_dlarfb.o ../org2r/src/my_dorgkr.o ../org2r/src/lumm.o ../org2r/src/dtvt.o ../org2r/src/dtrmmoop.o ../org2r/src/my_dlarft.o 
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

test_v7.exe: test_my_dorgqr.o my_dorgqr_v7.o ../org2r/src/my_dorgkr.o ../org2r/src/lumm.o ../org2r/src/dtvt.o ../org2r/src/dtrmmoop.o ../org2r/src/my_dlarft.o 
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

test.exe: my_dorgqr.o test_my_dorgqr.o
	$(FC) $(CFLAGS) $^ $(LAPACKELIB) $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

timeDorgqr.exe: timeDorgqr.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver: driver.c
	$(CC) $(CFLAGS) $^ -o $@

# Object file creations

my_dlarfb.o: my_dlarfb.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_my_dorgqr.o: test.c
	$(CC) -c $(CFLAGS) -DUSE_MY_DORGQR $^ -o $@

timeSVD.o: timeSVD.c
	$(CC) -c $(CFLAGS) -DUSE_AMD=1 $^ -o $@

my_dgesvd.o: my_dgesvd.f
	$(FC) -c $(FFLAGS) $^ -o $@

timeDorgqr.o: timeDorgqr.c
	$(CC) -c $(CFLAGS) $^ -o $@

dgesvd.o: dgesvd.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr.o: my_dorgqr.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v0.o: my_dorgqr_v0.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v1.o: my_dorgqr_v1.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v2.o: my_dorgqr_v2.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v3.o: my_dorgqr_v3.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v4.o: my_dorgqr_v4.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v5.o: my_dorgqr_v5.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgqr_v7.o: my_dorgqr_v7.f
	$(FC) -c $(FFLAGS) $^ -o $@

test.o: test.c 
	$(CC) -c $(CFLAGS) -I$(LAPACKEINC) -I$(CBLASINC) $^ -o $@

dorgqr.o: dorgqr.f
	$(FC) -c $(CFLAGS) $^ -o $@
