pThardmin_array=(10 20 30 40)
JetR_array=(0.2 0.3 0.4 0.5 0.6)

for(( i=0; i <=3;  i++ ))
  do
    for(( j=0; j <=4; j++ ))
      do  
        root -b -l -q "maindriver_newton.C( 100000000, 39034, 100, ${JetR_array[$j]}, 0 , 1000 , 0 , ${pThardmin_array[$i]} , 1 , 0., kTRUE, kFALSE , kFALSE )"
      done
  done


