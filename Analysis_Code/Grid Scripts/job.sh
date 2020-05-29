#!/bin/bash 

#source /nics/a/proj/UTK0019/alicesw/alice-env.sh -n 1
#export PYTHIA6428=/nics/a/proj/UTK0019/alicesw/pythia/pythia6428
#export outdir=/nics/b/home/chughe26/scratch/jetbin1/$PBS_JOBID
export outdir=/lustre/haven/user/chughe26/scratch/jetbin1/$PBS_JOBID
mkdir $outdir

#cd $PBS_O_WORKDIR
#cd /nics/d/home/cnattras/scratch/B2

echo $PBS_TASKNUM
JOBID=`echo $PBS_JOBID | cut -d . -f 1`
#To have completely unique seeds, need to multiply by a number which is greater than the number of jobs.  Protecting in case multiple jobs start at once.
NEWTASKNUM=$((PBS_TASKNUM * 20000))
JOB_NUMBER=$((NEWTASKNUM + $JOBID))

echo $JOB_NUMBER 
pwd
mkdir $outdir/$JOB_NUMBER
cp /nics/b/home/chughe26/ML_Studies/ACF_Analysis_Code/BkgrLoad.h $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/ML_Studies/ACF_Analysis_Code/BkgrLoad.cxx $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/ML_Studies/ACF_Analysis_Code/maindriver_newton.C $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/ML_Studies/ACF_Analysis_Code/Toy_Model_ML_Study.C $outdir/$JOB_NUMBER/ 

echo $1
echo $PBS_TASKNUM
cp /lustre/haven/user/chughe26/scratch/BKGD_EVENTS/0_5_percent/v1_v5/FINAL_FILES/$1/$PBS_TASKNUM/*.root $outdir/$JOB_NUMBER/

cd $outdir/$JOB_NUMBER
#root -b -l -q  "maindriver_newton.C(100 , $JOB_NUMBER  , 100 , 0.2 , 0 , 1000 , 0  , 1 )"

#root -b -l -q "maindriver_newton.C( 100000000, $JOB_NUMBER, 100, 0.4, 0 , 1000 , 0 , 1, 0.,  kTRUE, kFALSE, kTRUE )"

pThardmin_array=(10 , 20 , 30 , 40)
JetR_array=(0.2 , 0.3, 0.4, 0.5, 0.6)

for(( i=0 , i <=3 , i++ ))
  do
    for(( j=0, j <=4 , j++ ))
      do  
        root -b -l -q "maindriver_newton.C( 100000000, $JOB_NUMBER, 100, ${JetR_array[$j]}, 0 , 1000 , 0 , ${pThardmin_array[$i]} , 1 , 0., kTRUE, kFALSE , kTRUE )"
      done
  done
# num_pythia events / seed / pythia tune / jet radius / Harmonic Flag / DCA / centrality bin / pT hard min / Number of BKGD Files (BKGD events / 100 ) / constit cut / Data Flag (legacy to be removed) / Runtime Stats Flag / Run on Grid Flag
#rm -r $outdir/$JOB_NUMBER

