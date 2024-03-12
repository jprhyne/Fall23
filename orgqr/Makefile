include make.inc

2rObjs := $(patsubst %.f,%.o,$(wildcard ../org2r/src/*.f))
objects := $(patsubst %.f,%.o,$(wildcard ./src/*.f))
# Group definitions
default: srcFiles test.exe driver timeSVD.exe

srcFiles:
	+$(MAKE) -C src

clean: cleanSrc
	rm driver *exe

cleanSrc:
	+$(MAKE) -C src clean

# Executable compilations
timeSVD_opt.exe: timeSVD.o src/my_dorgqr_v1.o src/my_dgesvd.o
	$(FC) $(CFLAGS) $^ $(LAPACKELIB) $(OPTBLAS) -o $@

timeSVD.exe: timeSVD.o src/my_dorgqr.o src/my_dgesvd.o referenceFiles/dgesvd.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) $(LAPACKELIB) -o $@

test.exe: test.o $(2rObjs) $(objects)
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

timeDorgqr.exe: timeDorgqr.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@

driver: driver.c
	$(CC) $(CFLAGS) $^ -o $@

# Object file creations

test.o: test.c
	$(CC) -c $(CFLAGS) $^ -o $@

timeSVD.o: timeSVD.c
	$(CC) -c $(CFLAGS) -DUSE_AMD=1 $^ -o $@

timeDorgqr.o: timeDorgqr.c
	$(CC) -c $(CFLAGS) $^ -o $@
