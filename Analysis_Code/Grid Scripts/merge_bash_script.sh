#!/bin/bash
#for (( i=1949417; i<=1959470; i++ ))
#for (( i=1896303; i<=1896304; i++ ))
for (( i=3628336; i<=3628362; i++ ))
do
   JOB_NUM=$((i))
   DIR_STR=$"/lustre/haven/user/chughe26/scratch/jetbin1/$JOB_NUM"".apollo-acf"
   echo $DIR_STR
   FILE_STR11=$"$JOB_NUM""-jet-=0.2-pThard=10.root"
   CSV_STR11=$"$JOB_NUM""-jet-=0.2-pThard=10.csv"
   FILE_STR12=$"$JOB_NUM""-jet-=0.2-pThard=20.root"
   CSV_STR12=$"$JOB_NUM""-jet-=0.2-pThard=20.csv"
   FILE_STR13=$"$JOB_NUM""-jet-=0.2-pThard=30.root"
   CSV_STR13=$"$JOB_NUM""-jet-=0.2-pThard=30.csv"
   FILE_STR14=$"$JOB_NUM""-jet-=0.2-pThard=40.root"
   CSV_STR14=$"$JOB_NUM""-jet-=0.2-pThard=40.csv"
   FILE_STR21=$"$JOB_NUM""-jet-=0.3-pThard=10.root"
   CSV_STR21=$"$JOB_NUM""-jet-=0.3-pThard=10.csv"
   FILE_STR22=$"$JOB_NUM""-jet-=0.3-pThard=20.root"
   CSV_STR22=$"$JOB_NUM""-jet-=0.3-pThard=20.csv"
   FILE_STR23=$"$JOB_NUM""-jet-=0.3-pThard=30.root"
   CSV_STR23=$"$JOB_NUM""-jet-=0.3-pThard=30.csv"
   FILE_STR24=$"$JOB_NUM""-jet-=0.3-pThard=40.root"
   CSV_STR24=$"$JOB_NUM""-jet-=0.3-pThard=40.csv"
   FILE_STR31=$"$JOB_NUM""-jet-=0.4-pThard=10.root"
   CSV_STR31=$"$JOB_NUM""-jet-=0.4-pThard=10.csv"
   FILE_STR32=$"$JOB_NUM""-jet-=0.4-pThard=20.root"
   CSV_STR32=$"$JOB_NUM""-jet-=0.4-pThard=20.csv"
   FILE_STR33=$"$JOB_NUM""-jet-=0.4-pThard=30.root"
   CSV_STR33=$"$JOB_NUM""-jet-=0.4-pThard=30.csv"
   FILE_STR34=$"$JOB_NUM""-jet-=0.4-pThard=40.root"
   CSV_STR34=$"$JOB_NUM""-jet-=0.4-pThard=40.csv"
   FILE_STR41=$"$JOB_NUM""-jet-=0.5-pThard=10.root"
   CSV_STR41=$"$JOB_NUM""-jet-=0.5-pThard=10.csv"
   FILE_STR42=$"$JOB_NUM""-jet-=0.5-pThard=20.root"
   CSV_STR42=$"$JOB_NUM""-jet-=0.5-pThard=20.csv"
   FILE_STR43=$"$JOB_NUM""-jet-=0.5-pThard=30.root"
   CSV_STR43=$"$JOB_NUM""-jet-=0.5-pThard=30.csv"
   FILE_STR44=$"$JOB_NUM""-jet-=0.5-pThard=40.root"
   CSV_STR44=$"$JOB_NUM""-jet-=0.5-pThard=40.csv"
   FILE_STR51=$"$JOB_NUM""-jet-=0.6-pThard=10.root"
   CSV_STR51=$"$JOB_NUM""-jet-=0.6-pThard=10.csv"
   FILE_STR52=$"$JOB_NUM""-jet-=0.6-pThard=20.root"
   CSV_STR52=$"$JOB_NUM""-jet-=0.6-pThard=20.csv"
   FILE_STR53=$"$JOB_NUM""-jet-=0.6-pThard=30.root"
   CSV_STR53=$"$JOB_NUM""-jet-=0.6-pThard=30.csv"
   FILE_STR54=$"$JOB_NUM""-jet-=0.6-pThard=40.root"
   CSV_STR54=$"$JOB_NUM""-jet-=0.6-pThard=40.csv"
if [ -d $DIR_STR ]
then
hadd $FILE_STR11 $DIR_STR/3*/*0.2*10*.root
hadd $FILE_STR12 $DIR_STR/3*/*0.2*20*.root
hadd $FILE_STR13 $DIR_STR/3*/*0.2*30*.root
hadd $FILE_STR14 $DIR_STR/3*/*0.2*40*.root
hadd $FILE_STR21 $DIR_STR/3*/*0.3*10*.root
hadd $FILE_STR22 $DIR_STR/3*/*0.3*20*.root
hadd $FILE_STR23 $DIR_STR/3*/*0.3*30*.root
hadd $FILE_STR24 $DIR_STR/3*/*0.3*40*.root
hadd $FILE_STR31 $DIR_STR/3*/*0.4*10*.root
hadd $FILE_STR32 $DIR_STR/3*/*0.4*20*.root
hadd $FILE_STR33 $DIR_STR/3*/*0.4*30*.root
hadd $FILE_STR34 $DIR_STR/3*/*0.4*40*.root
hadd $FILE_STR41 $DIR_STR/3*/*0.5*10*.root
hadd $FILE_STR42 $DIR_STR/3*/*0.5*20*.root
hadd $FILE_STR43 $DIR_STR/3*/*0.5*30*.root
hadd $FILE_STR44 $DIR_STR/3*/*0.5*40*.root
hadd $FILE_STR51 $DIR_STR/3*/*0.6*10*.root
hadd $FILE_STR52 $DIR_STR/3*/*0.6*20*.root
hadd $FILE_STR53 $DIR_STR/3*/*0.6*30*.root
hadd $FILE_STR54 $DIR_STR/3*/*0.6*40*.root
cat $DIR_STR/3*/*0.2*10*.csv > $CSV_STR11
cat $DIR_STR/3*/*0.2*20*.csv > $CSV_STR12
cat $DIR_STR/3*/*0.2*30*.csv > $CSV_STR13
cat $DIR_STR/3*/*0.2*40*.csv > $CSV_STR14	
cat $DIR_STR/3*/*0.3*10*.csv > $CSV_STR21
cat $DIR_STR/3*/*0.3*20*.csv > $CSV_STR22
cat $DIR_STR/3*/*0.3*30*.csv > $CSV_STR23
cat $DIR_STR/3*/*0.3*40*.csv > $CSV_STR24
cat $DIR_STR/3*/*0.4*10*.csv > $CSV_STR31
cat $DIR_STR/3*/*0.4*20*.csv > $CSV_STR32
cat $DIR_STR/3*/*0.4*30*.csv > $CSV_STR33
cat $DIR_STR/3*/*0.4*40*.csv > $CSV_STR34	
cat $DIR_STR/3*/*0.5*10*.csv > $CSV_STR41
cat $DIR_STR/3*/*0.5*20*.csv > $CSV_STR42
cat $DIR_STR/3*/*0.5*30*.csv > $CSV_STR43
cat $DIR_STR/3*/*0.5*40*.csv > $CSV_STR44
cat $DIR_STR/3*/*0.6*10*.csv > $CSV_STR51
cat $DIR_STR/3*/*0.6*20*.csv > $CSV_STR52
cat $DIR_STR/3*/*0.6*30*.csv > $CSV_STR53
cat $DIR_STR/3*/*0.6*40*.csv > $CSV_STR54
fi
done

hadd merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-10.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.2*10*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.2*10.csv > merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-10.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-20.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.2*20*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.2*20.csv > merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-20.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-30.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.2*30*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.2*30.csv > merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-30.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-40.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.2*40*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.2*40.csv > merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-40.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.3-pThardmin-10.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.3*10*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.3*10.csv > merged-ML-output-LOWSTATS-Rparam-0.3-pThardmin-10.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.3-pThardmin-20.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.3*20*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.3*20.csv > merged-ML-output-LOWSTATS-Rparam-0.3-pThardmin-20.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.3-pThardmin-30.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.3*30*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.3*30.csv > merged-ML-output-LOWSTATS-Rparam-0.3-pThardmin-30.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.3-pThardmin-40.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.3*40*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.3*40.csv > merged-ML-output-LOWSTATS-Rparam-0.3-pThardmin-40.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.4-pThardmin-10.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.4*10*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.4*10.csv > merged-ML-output-LOWSTATS-Rparam-0.4-pThardmin-10.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.4-pThardmin-20.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.4*20*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.4*20.csv > merged-ML-output-LOWSTATS-Rparam-0.4-pThardmin-20.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.4-pThardmin-30.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.4*30*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.4*30.csv > merged-ML-output-LOWSTATS-Rparam-0.4-pThardmin-30.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.4-pThardmin-40.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.4*40*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.4*40.csv > merged-ML-output-LOWSTATS-Rparam-0.4-pThardmin-40.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.5-pThardmin-10.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.5*10*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.5*10.csv > merged-ML-output-LOWSTATS-Rparam-0.5-pThardmin-10.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.5-pThardmin-20.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.5*20*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.5*20.csv > merged-ML-output-LOWSTATS-Rparam-0.5-pThardmin-20.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.5-pThardmin-30.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.5*30*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.5*30.csv > merged-ML-output-LOWSTATS-Rparam-0.5-pThardmin-30.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.5-pThardmin-40.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.5*40*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.5*40.csv > merged-ML-output-LOWSTATS-Rparam-0.5-pThardmin-40.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.6-pThardmin-10.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.6*10*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.6*10.csv > merged-ML-output-LOWSTATS-Rparam-0.6-pThardmin-10.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.6-pThardmin-20.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.6*20*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.6*20.csv > merged-ML-output-LOWSTATS-Rparam-0.6-pThardmin-20.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.6-pThardmin-30.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.6*30*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.6*30.csv > merged-ML-output-LOWSTATS-Rparam-0.6-pThardmin-30.csv

hadd merged-ML-output-LOWSTATS-Rparam-0.6-pThardmin-40.root /lustre/haven/user/chughe26/scratch/jetbin1/*0.6*40*.root
cat /lustre/haven/user/chughe26/scratch/jetbin1/*0.6*40.csv > merged-ML-output-LOWSTATS-Rparam-0.6-pThardmin-40.csv
