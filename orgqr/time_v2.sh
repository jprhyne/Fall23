#!/bin/env bash
#This file will time driver using test.exe as well as test_optBlas.exe for the reference and optimized BLAS/LAPACK routines resp.

# First, we run make clean to clear out the directory
#make clean

# Next, we run make to compile our 
#make


# Now, we will run test.exe for varying m,n,k files only grabbing the timing
echo "Tall and skinny matrices"
for (( m=100000; m<=1000000; m+=200000 ))
do
    for (( n=10; n<=400; n+=10 ))
    do
        echo "m=$m n=$n k=$n"
        echo "my_dorgqr ran 15 times"
        for (( l=1; l<=15; l+=1 ))
        do
            ./test_v2.exe -m $m -n $n -k $n -t
        done
    done
done
echo "Square matrices"
for (( m=500; m<=5000; m+=500 ))
do
    n=$m
    for (( k=500; k<=n; k+=500 ))
    do
        echo "m=$m n=$n k=$k"
        echo "my_dorgqr ran 15 times"
        for (( l=1; l<=15; l+=1 ))
        do
            ./test_v2.exe -m $m -n $n -k $k -t
        done 
    done
done
