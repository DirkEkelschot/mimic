for i in `seq $1 $2`;
    do
    touch post.scr;
    echo "pob_init $i";
    echo "pob_init $i" >> post.scr;
    echo "pob_write_fvs file=fvs_"$i".vtk">> post.scr;
#   echo "pob_write_fvs file=fvs_"$i".dat">> post.scr;
    echo "quit">> post.scr;
    us3d-postpar -D --script=post.scr
    rm post.scr
#   pplt fvs_$i.dat
done
