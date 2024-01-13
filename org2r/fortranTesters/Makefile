include ../make.inc

# Groups
all: objectFiles driver_dtvt.exe driver_trmm.exe driver_dorgkr.exe driver_lumm.exe driver_dlarft.exe

clean:
	rm *.o *.exe

objectFiles: 
	$(MAKE) -C ../src all

driver_dtvt.exe: driver_dtvt.o test_dtvt.o  ../src/dtrmmoop.o ../src/dtvt.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

driver_trmm.exe: driver_trmm.o testdtrmmoop.o ../src/dtrmmoop.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

driver_dorgkr.exe: driver_dorgkr.o test_dorgkr.o ../src/my_dorgkr.o ../src/lumm.o ../src/dtrmmoop.o ../src/dtvt.o  ../src/my_dlarft.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

driver_dlarft.exe: driver_dlarft.o test_dlarft.o ../src/my_dorgkr.o ../src/lumm.o ../src/dtrmmoop.o ../src/dtvt.o ../src/my_dlarft.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

driver_lumm.exe:   driver_lumm.o test_lumm.o ../src/lumm.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) -o $@

driver_dtvt.o: driver_dtvt.c
	$(CC) -c $(CFLAGS) $^ -o $@

driver_dlarft.o: driver_dlarft.c
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

test_dlarft.o: test_dlarft.f
	$(FC) -c $(FFLAGS) $^ -o $@

test_lumm.o: test_lumm.f
	$(FC) -c $(FFLAGS) $^ -o $@

testdtrmmoop.o: testdtrmmoop.f
	$(FC) -c $(FFLAGS) $^ -o $@

