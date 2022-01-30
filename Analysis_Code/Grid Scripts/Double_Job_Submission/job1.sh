#!/bin/bash 

##create calculation sub-directories######################################################################
#"To have completely unique seeds, need to multiply by a number which is greater than the number of jobs."
#"Protecting in case multiple jobs start at once."
export outdir=${PBS_O_WORKDIR}/$PBS_JOBID
#export outdir=/lustre/haven/user/chughe26/scratch/jetbin6/$PBS_JOBID
mkdir $outdir
echo $PBS_TASKNUM
JOBID=`echo $PBS_JOBID | cut -d . -f 1`
NEWTASKNUM=$((PBS_TASKNUM * 20000))
JOB_NUMBER=$((NEWTASKNUM + $JOBID))
echo $JOB_NUMBER 
pwd
mkdir $outdir/$JOB_NUMBER

##copy input codes#############################################################################
cp /nics/b/home/chughe26/ML_Studies/ACF_Analysis_Code/BkgrLoad.h $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/ML_Studies/ACF_Analysis_Code/BkgrLoad.cxx $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/ML_Studies/ACF_Analysis_Code/maindriver_newton.C $outdir/$JOB_NUMBER/
cp /nics/b/home/chughe26/ML_Studies/ACF_Analysis_Code/Toy_Model_ML_Study.C $outdir/$JOB_NUMBER/ 

##copy input files (HI ONLY !!!)###############################################################
echo $1
echo $PBS_TASKNUM
cp /lustre/haven/user/chughe26/scratch/BKGD_EVENTS/0_5_percent/v1_v5/FINAL_FILES/$1/$PBS_TASKNUM/*.root $outdir/$JOB_NUMBER/

##load modules#################################################################################
## Do NOT change the order of the following files
cd $outdir/$JOB_NUMBER
export PATH=$PATH:/lustre/haven/proj/UTK0019/python-virtual-environments/env/bin
export ALIBUILD_WORK_DIR=/lustre/haven/proj/UTK0019/aliroot5/alice/sw
eval "`alienv shell-helper`"
#alienv load AliPhysics/latest
alienv --no-refresh load AliPhysics/latest


##initiate calculation#########################################################################
#root -b -l -q  "maindriver_newton.C(100 , $JOB_NUMBER  , 100 , 0.2 , 0 , 1000 , 0  , 1 )"

#root -b -l -q "maindriver_newton.C( 100000000, $JOB_NUMBER, 100, 0.4, 0 , 1000 , 0 , 1, 0.,  kTRUE, kFALSE, kTRUE )"

pThardmin_array=(10.0 20.0 30.0)
JetR_array=(0.2 0.3 0.4 0.5 0.6)

for(( i=0; i <=2; i++ ))
  do
    for(( j=0; j <=4; j++ ))
      do  
        echo '________________________________________________________________________________________________'
        var=$(printf '\npThardmin = %s, Rparam = %s\n' "${pThardmin_array[$i]}" "${JetR_array[$j]}")
        echo $var
        echo '________________________________________________________________________________________________'
        root -b -l -q "maindriver_newton.C( 100000000, $JOB_NUMBER, 356, ${JetR_array[$j]}, 0.1 , 0 , 1000 , 0 , ${pThardmin_array[$i]} , 1 , 0., kTRUE, kFALSE , kTRUE )"
      done
  done
# num_pythia events / seed / pythia tune / jet radius / subjet radius / Harmonic Flag / DCA / centrality bin / pT hard min / Number of BKGD Files (BKGD events / 100 ) / constit cut / Data Flag (legacy to be removed) / Runtime Stats Flag / Run on Grid Flag

################################################################################################
#rm -r $outdir/$JOB_NUMBER

