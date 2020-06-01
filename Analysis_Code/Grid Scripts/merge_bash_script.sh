#!/bin/bash
#for (( i=1949417; i<=1959470; i++ ))
#for (( i=1896303; i<=1896304; i++ ))

pThardmin_array=(10.0 20.0 30.0 40.0)
JetR_array=(0.2 0.3 0.4 0.5 0.6)

for (( i=3657797; i<=3657823; i++ ))
do
   JOB_NUM=$((i))
   DIR_STR=$"/lustre/haven/user/chughe26/scratch/jetbin1/$JOB_NUM"".apollo-acf"
   echo $DIR_STR


if [ -d $DIR_STR ]
then

for(( w=0; w <=3; w++ ))
  do
    for(( x=0; x <=4; x++ ))
      do  
     
        ROOTFILE_STR=$(printf '%s-Rparam-=-%s-pThardmin-=-%s-.root' "$JOB_NUM" "${JetR_array[$x]}" "${pThardmin_array[$w]}")
        CSVFILE_STR=$(printf '%s-Rparam-=-%s-pThardmin-=-%s-.csv' "$JOB_NUM" "${JetR_array[$x]}" "${pThardmin_array[$w]}")

        hadd $ROOTFILE_STR $DIR_STR/3*/*${JetR_array[$x]}*${pThardmin_array[$w]}*.root
        cat $DIR_STR/3*/*${JetR_array[$x]}*${pThardmin_array[$w]}*.csv > $CSVFILE_STR
      done
  done

fi
done


for(( r=0; r <=3; r++ ))
  do
    for(( s=0; s <=4; s++ ))
      do  

        MERGE_ROOTFILE_STR=$(printf 'merged-ML-output-LOWSTATS-Rparam-%s-pThardmin-%s.root' "${JetR_array[$s]}" "${pThardmin_array[$r]}")
        MERGE_CSVFILE_STR=$(printf 'merged-ML-output-LOWSTATS-Rparam-%s-pThardmin-%s.csv' "${JetR_array[$s]}" "${pThardmin_array[$r]}")       

        hadd $MERGE_ROOTFILE_STR /lustre/haven/user/chughe26/scratch/jetbin1/3*${JetR_array[$s]}*${pThardmin_array[$r]}*.root
        cat $DIR_STR/3*${JetR_array[$s]}*${pThardmin_array[$r]}*.csv > $MERGE_CSVFILE_STR
      done
  done

