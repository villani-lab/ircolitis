#!/usr/bin/env bash

dirs=(
a12_4_4_b5_1_3 
a12_4_4_m3_2  
a12_4_4_t4_cd4_2_2
a12_4_4_t4_cd8_1_2
blood2_bcell5  
blood2_myeloid5  
blood2_tcell5_cd4_5
blood2_tcell5_cd8_5  
n3_2
)


for d in ${dirs[*]}
do
  command ln -s ../a20/${d}_html/${d} .
done
