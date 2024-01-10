include make.inc
# Group definitions
default: overallTest.exe overallTest_optBlas.exe test_dtrmmoop.exe test_dtrmmoop_optBlas.exe driver_trmm.exe driver_dorgkr.exe driver_lumm.exe driver_dtvt.exe

clean:
	rm *o *exe

# Executable compilations
overallTest_optBlas.exe: overallTest_opt.o
	$(FC) -D USE_AOCL $(CFLAGS) $^ $(OPTBLAS) -o $@

overallTest.exe: overallTest.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

test_dtrmmoop.exe: test_dtrmmoop.o  dtrmmoop.o
	$(FC) $^ $(LAPACKLIB) $(BLASLIB) -o $@

test_dtrmmoop_optBlas.exe: test_dtrmmoop_opt.o dtrmmoop.o
	$(FC) $^ $(OPTBLAS) -o $@

driver_dtvt.exe: driver_dtvt.o test_dtvt.o  dtrmmoop.o dtvt.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

driver_trmm.exe: driver_trmm.o testdtrmmoop.o dtrmmoop.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

driver_dorgkr.exe: driver_dorgkr.o test_dorgkr.o my_dorgkr.o lumm.o dtrmmoop.o dtvt.o 
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

driver_lumm.exe:   driver_lumm.o test_lumm.o lumm.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

# Object file creations

driver_dtvt.o: driver_dtvt.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_trmm.o: driver_trmm.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_lumm.o: driver_lumm.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_dorgkr.o: driver_dorgkr.c
	$(CC) -c $(CFLAGS) $^ -o $@

test_dtvt.o: test_dtvt.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_dorgkr.o: test_dorgkr.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_lumm.o: test_lumm.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dorgkr.o: my_dorgkr.f
	$(FC) -c $(FFLAGS) $^ -o $@

lumm.o: lumm.f
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
