include make.inc
# Group definitions
default: overallTest.exe overallTest_optBlas.exe test_dtrmmoop.exe test_dtrmmoop_optBlas.exe driver.exe

clean:
	rm *o *exe

# Executable compilations
overallTest_optBlas.exe: overallTest_opt.o overallTest.c
	$(FC) -D USE_AOCL $(CFLAGS) overallTest_opt.o $(OPTBLAS) -o $@

overallTest.exe: overallTest.o overallTest.c
	$(FC) $(CFLAGS) overallTest.o $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

test_dtrmmoop.exe: test_dtrmmoop.o test_dtrmmoop.c dtrmmoop.o
	$(FC) test_dtrmmoop.o dtrmmoop.o $(LAPACKLIB) $(BLASLIB) -o $@

test_dtrmmoop_optBlas.exe: test_dtrmmoop_opt.o test_dtrmmoop.c dtrmmoop.o
	$(FC) test_dtrmmoop_opt.o dtrmmoop.o $(OPTBLAS) -o $@

test_dtvt.exe: test_dtvt.o dtvt.o test_dtvt.c
	$(FC) $(CFLAGS) test_dtvt.o dtvt.o -o $@

driver.exe: driver.o testdtrmmoop.o dtrmmoop.o
	$(FC) $(CFLAGS) driver.o testdtrmmoop.o dtrmmoop.o $(LAPACKLIB) $(BLASLIB) -o $@

# Object file creations

driver.o: driver.c
	$(CC) -c $(CFLAGS) $^ -o $@

testdtrmmoop.o: testdtrmmoop.f
	$(FC) -c $(FFLAGS) $^ -o $@

overallTest.o: overallTest.c 
	$(CC) -c $(CFLAGS) $^ -o $@

overallTest_opt.o: overallTest.c 
	$(CC) -D USE_AOCL -c $(CFLAGS) -I$(OPTBLASINC) -I$(OPTLAPACKINC) $^ -o $@

test_dtrmmoop.o: test_dtrmmoop.c
	$(FC) -c $(CFLAGS) $^ -o $@

test_dtrmmoop_opt.o: test_dtrmmoop.c
	$(CC) -D USE_AOCL -c $(CFLAGS) $^ -o $@

dtrmmoop.o: dtrmmoop.f
	$(FC) -c $(FFLAGS) $^ -o $@

dtvt.o: dtvt.f
	$(FC) -c $(FFLAGS) $^ -o $@
