#!/bin/bash

echo "______________________________________"
echo "STARTING pT_hard_min_10_20_30"

for (( i=2474; i<=2500; i++ ))
do 
   echo $i 
   ./qsubmit1.sh $i
done

echo "______________________________________"
echo "STARTING pT_hard_min_40_60_80"

for (( i=2474; i<=2500; i++ ))
do
   echo $i 
   ./qsubmit2.sh $i
done



#valid numbers 2474 - 2500 OR 7540 - 23704
