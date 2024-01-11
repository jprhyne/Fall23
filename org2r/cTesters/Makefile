include ../make.inc
# Groups

all: objectFiles overallTest_optBlas.exe overallTest.exe test_dtrmmoop.exe test_dtrmmoop_optBlas.exe

clean:
	rm *.o *.exe

objectFiles:
	$(MAKE) -C ../src all

overallTest_optBlas.exe: overallTest_opt.o ../src/my_dorgkr.o ../src/lumm.o ../src/dtvt.o ../src/dtrmmoop.o
	$(FC) -D USE_AOCL $(CFLAGS) $^ $(OPTBLAS) -o $@

overallTest.exe: overallTest.o ../src/my_dorgkr.o ../src/lumm.o ../src/dtvt.o ../src/dtrmmoop.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(CBLASLIB) $(BLASLIB) -o $@

test_dtrmmoop.exe: test_dtrmmoop.o  ../src/dtrmmoop.o
	$(FC) $^ $(LAPACKLIB) $(BLASLIB) -o $@

test_dtrmmoop_optBlas.exe: test_dtrmmoop_opt.o ../src/dtrmmoop.o
	$(FC) $^ $(OPTBLAS) -o $@

overallTest.o: overallTest.c 
	$(CC) -c $(CFLAGS) $^ -o $@

overallTest_opt.o: overallTest.c 
	$(CC) -D USE_AOCL -c $(CFLAGS) -I$(OPTBLASINC) -I$(OPTLAPACKINC) $^ -o $@

test_dtrmmoop.o: test_dtrmmoop.c
	$(FC) -c $(CFLAGS) $^ -o $@

test_dtrmmoop_opt.o: test_dtrmmoop.c
	$(CC) -D USE_AOCL -c $(CFLAGS) $^ -o $@

