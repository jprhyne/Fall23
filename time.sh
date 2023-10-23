#!/bin/env bash
#This file will time driver using test.exe as well as test_optBlas.exe for the reference and optimized BLAS/LAPACK routines resp.

# First, we run make clean to clear out the directory
#make clean

# Next, we run make to compile our 
#make

# Now, we will run test.exe for varying m,n,k files only grabbing the GFlop/s 
for (( m=50; m<=850; m+=100 ))
do
    for (( n=50; n<=m; n+=100 ))
    do
        for (( k=10; k<=n; k+=100 ))
        do
            echo "m=$m n=$n k=$k"
            echo "Reference dorgqr"
            ./test.exe -m $m -n $n -k $k -t

            echo "Optimized dorgqr"
            echo "my_dorgqr"
            ./test_optBlas.exe -m $m -n $n -k $k -t
            sleep 1
        done
    done
done
