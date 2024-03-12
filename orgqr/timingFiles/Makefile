include make.inc

2rObjs := $(patsubst %.f,%.o,$(wildcard ../../org2r/src/*.f))
objects := $(patsubst %.f,%.o,$(wildcard ../src/*.f))

default: srcFiles timeSVD.exe timeSVD_opt.exe timeDorgqr.exe timeDlarfb.exe

srcFiles:
	+$(MAKE) -C ../src

# Executable compilations
timeSVD_opt.exe: timeSVD.o src/my_dorgqr_v1.o src/my_dgesvd.o
	$(FC) $(CFLAGS) $^ $(LAPACKELIB) $(OPTBLAS) -o $@

timeSVD.exe: timeSVD.o src/my_dorgqr.o src/my_dgesvd.o referenceFiles/dgesvd.o
	$(FC) $(CFLAGS) $^ $(LAPACKLIB) $(BLASLIB) $(LAPACKELIB) -o $@

timeDorgqr.exe: timeDorgqr.o
	$(FC) $(CFLAGS) $^ $(OPTBLAS) -o $@
# Object file creations
timeSVD.o: timeSVD.c
	$(CC) -c $(CFLAGS) -DUSE_AMD=1 $^ -o $@

timeDorgqr.o: timeDorgqr.c
	$(CC) -c $(CFLAGS) $^ -o $@

timeDlarfb.o: timeDlarfb.c
	$(CC) -c $(CFLAGS) $^ -o $@
