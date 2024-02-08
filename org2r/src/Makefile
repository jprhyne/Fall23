include ../make.inc

# Groups
all: my_dlarft.o my_dorgkr.o lumm.o dtrmmoop.o dtvt.o

clean:
	rm *.o

my_dorgkr.o: my_dorgkr.f
	$(FC) -c $(FFLAGS) $^ -o $@

my_dlarft.o: my_dlarft.f
	$(FC) -c $(FFLAGS) $^ -o $@

lumm.o: lumm.f
	$(FC) -c $(FFLAGS) $^ -o $@

dtrmmoop.o: dtrmmoop.f
	$(FC) -c $(FFLAGS) $^ -o $@

dtvt.o: dtvt.f
	$(FC) -c $(FFLAGS) $^ -o $@

dorg2r.o: dorg2r.f
	$(FC) -c $(FFLAGS) $^ -o $@
