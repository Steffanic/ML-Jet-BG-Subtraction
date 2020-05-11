#!/bin/bash
#change indexes to go over the job you want to merge upon finishing
for (( i=2715375; i<=2715882; i++ ))
do
   JOB_NUM=$((i))
   DIR_STR=$"/lustre/haven/user/chughe26/scratch/jetbin1/$JOB_NUM"".apollo-acf"
   echo $DIR_STR
   FILE_STR=$"$JOB_NUM"".root"
   CSV_STR=$"$JOB_NUM"".csv"
if [ -d $DIR_STR ]
then
hadd $FILE_STR $DIR_STR/3*/*.root
cat $DIR_STR/3*/*.csv > $CSV_STR
fi
done
hadd merged_ML_output_LOW_STATS.root /lustre/haven/user/chughe26/scratch/jetbin1/3*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*.csv > merged_ML_output_LOW_STATS.csv

