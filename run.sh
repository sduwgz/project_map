run map
./map -m data/input2.bnx -g data/ref.cmp 1> data/debugdata/map_result 2> data/debugdata/map_pos

compare
cat map_result | awk '{print $1,$4,$5}' | python py/script/compare.py huada_map > cpmpare_new_out

python py/script/transform.py map_result | python py/script/compare.py huada_map

cat result/map_result_01 | awk '{if(!i) {s=$0;ss=$2;}else if($2>ss) print $0;else print s;i=!i}' | awk '{if($2<12) print $0}' > result/map_result_low_score
cat result/map_result_2600 | python py/script/compare_tool3.py result/map_result_low_score | python py/script/compare.py huada_map
cat result/map_result_01 | python py/script/compare.py huada_map | python py/script/compare_tool4.py result/map_pos_01 | awk '{if($6-$5>-2000 && $6-$5<2000) print $5-$6}' > temp
cat poisson_data | awk '{print int($1*100+0.5)}' | awk 'BEGIN{count=0}{if($1==1) count++;}END{print count}'

get the higher mapping between the mole and its reversed one.
cat result/standard_data_result | awk 'BEGIN{i=1}{if(i){s=$0;ss=$2;i=!i}else{if(ss>$2) print s;else print $0;i=!i}}' | python py/script/compare.py data/huada_map_6943
cat result/standard_data_result | awk '{i=!i;if(i){s=$0;ss=$2;}else{if(ss>$2){if(ss-$2<5)print $0; print s;}else{if($2-ss<5) print s; print $0;}}}' | python py/script/compare.py data/huada_map_6943
cat check_result/result_18446 |  python py/script/stairs.py data/ref_filter.cmap data/input2.bnx 18446

