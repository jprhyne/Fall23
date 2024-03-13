#!/bin/env bash

# Fix n and k
# n is a decently sized number of columns, and k is the blocksize. We are sticking with
# 32 in this as we are concerned with performance inside dorgqr for dlarfb
for (( m=50000; m<= 500000; m+=20000 ))
do 
    echo "m=$m"
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./timeDlarfb.exe -m $m -n 30000 -k 32
    done
done
