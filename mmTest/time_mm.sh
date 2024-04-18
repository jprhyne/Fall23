for (( m=1000; m<=10000; m+=1000 ));
do
    echo "m=$m"
    n=$m
    k=$m
    echo "Testing file ran 10 times"
    for (( l=1; l<=10; l+=1 ))
    do
        ./test.exe -t -m $m -n $n -k $k
    done
done
