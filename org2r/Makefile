include make.inc
# Group definitions
#default: overallTest.exe overallTest_optBlas.exe test_dtrmmoop.exe test_dtrmmoop_optBlas.exe driver_trmm.exe driver_dorgkr.exe driver_lumm.exe driver_dtvt.exe

all: 
	+$(MAKE) -C src all
	+$(MAKE) -C fortranTesters all
	+$(MAKE) -C cTesters all

clean:
	+$(MAKE) -C src clean
	+$(MAKE) -C fortranTesters clean
	+$(MAKE) -C cTesters clean

