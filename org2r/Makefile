include make.inc
# Group definitions
default: overallTest.exe overallTest_optBlas.exe test_dtrmmoop.exe test_dtrmmoop_optBlas.exe driver_trmm.exe driver_dorg2r.exe

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

driver_trmm.exe: driver_trmm.o testdtrmmoop.o dtrmmoop.o
	$(FC) $(CFLAGS) driver_trmm.o testdtrmmoop.o dtrmmoop.o $(LAPACKLIB) $(BLASLIB) -o $@

driver_dorg2r.exe: driver_dorg2r.o test_dorg2r.o my_dorg2r_cheat.o lumm_cheat.o dtvt_cheat.o dtrmmoop.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

# Object file creations

driver_trmm.o: driver_trmm.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_dorg2r.o: driver_dorg2r.c
	$(CC) -c $(CFLAGS) $^ -o $@

test_dorg2r.o: test_dorg2r.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorg2r_cheat.o: my_dorg2r_cheat.f
	$(CC) -c $(CFLAGS) $^ -o $@

lumm_cheat.o: lumm_cheat.f
	$(FC) -c $(FFLAGS) $^ -o $@

dtvt_cheat.o: dtvt_cheat.f
	$(FC) -c $(FFLAGS) $^ -o $@

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
