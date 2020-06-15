pThardmin_array=(10.0 20.0 30.0 40.0)
JetR_array=(0.2 0.3 0.4 0.5 0.6)

for(( i=0; i <=3;  i++ ))
  do
    for(( j=0; j <=4; j++ ))
      do  
        echo '________________________________________________________________________________________________'
        var=$(printf '\npThardmin = %s, Rparam = %s\n' "${pThardmin_array[$i]}" "${JetR_array[$j]}")
        echo $var
        echo '________________________________________________________________________________________________'
        root -b -l -q "maindriver_newton.C( 100000000, 39034, 356, ${JetR_array[$j]}, 0 , 1000 , 0 , ${pThardmin_array[$i]} , 1 , 0., kTRUE, kTRUE , kFALSE )"
      done
  done


