for (( m=500000; m<5000000; m+=200000 ))
do
    echo "m=$m"
    echo "Testing file ran 5 times"
    for (( l=1; l<= 5; l+=1 ))
    do
        ./timeDlarft.exe -t -m $m
    done
done

# fix m and vary n
m=500000
for (( n=32; n<3201; n+=32 ))
do
    echo "n=$n"
    echo "Teting file ran 5 times"
    for (( l=1; l<= 5; l+=1 ))
    do
        ./timeDlarft.exe -t -m $m -n $n
    done
done
