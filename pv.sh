for((i=0;i<1000;i++))
do
    python py/script/p_value.py 4600000 > data/p_test_ref/ref.cmap
    ./map -m data/p_test_mole -g data/p_test_ref/ref.cmap -o result/p_value/map_pos
done
