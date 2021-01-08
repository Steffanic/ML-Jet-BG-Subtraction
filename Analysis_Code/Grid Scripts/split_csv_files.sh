splitCsv() {
    HEADER=$(head -1 $1)
     echo $HEADER
    if [ -n "$2" ]; then
        CHUNK=$2
    else 
        CHUNK=1000
    fi

    TEMP=$1
    NAME=${TEMP%.csv}

    echo $NAME

    split -dl $CHUNK --additional-suffix .csv $1 $NAME"_MINI_"

    count=0
    for i in $NAME"_MINI_"*csv; do
    if [ $count -gt 0 ]; then
      echo -e "$HEADER\n$(cat $i)" > $i
    fi
    ((count++))
    done

}

pThardmin_array=(10.0 20.0 30.0 40.0)
JetR_array=(0.2 0.3 0.4 0.5 0.6)

MASTER_DIR_STR=$"/lustre/haven/user/chughe26/scratch/jetbin1/LOWSTATS/"
mkdir $MASTER_DIR_STR


for(( w=0; w <=3; w++ ))
  do
    for(( x=0; x <=4; x++ ))
      do  

        echo "Rparam = ${JetR_array[$x]}, pThardmin = ${pThardmin_array[$w]}"

        DIR_STR=$(printf 'Rparam-=-%s-pThardmin-=%s' "${JetR_array[$x]}" "${pThardmin_array[$w]}")

        mkdir $MASTER_DIR_STR/$DIR_STR

        cd $MASTER_DIR_STR/$DIR_STR
     
        MERGE_CSVFILE_STR=$(printf 'merged-ML-output-LOWSTATS-Rparam-%s-pThardmin-%s.csv' "${JetR_array[$x]}" "${pThardmin_array[$w]}")

        cp "/lustre/haven/user/chughe26/scratch/jetbin1/$MERGE_CSVFILE_STR" .
        splitCsv $MERGE_CSVFILE_STR 50000
        rm -rf $MERGE_CSVFILE_STR
      done
  done

cd /lustre/haven/user/chughe26/scratch/jetbin1
echo "CLEANING UP THE EXCESS/INTERMEDIATE FILES NOW !!!"
rm -rf 3*.csv
rm -rf 3*.root
echo "TAR-ING ALL THE REMANING ROOT FILES !"
tar -czvf ROOTFILES-LOWSTATS.tar.gz *.root
echo "TAR-ING ALL THE REMAINING CSV FILES !"
tar -czvf CSVFILES-LOWSTATS.tar.gz *.csv LOWSTATS/
