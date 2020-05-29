#!/bin/bash
for (( i=2474; i<=2500; i++ ))
do 
   echo $i 
   source submit.sh $i
done

#valid numbers 2474 - 2500 OR 7540 - 17538
