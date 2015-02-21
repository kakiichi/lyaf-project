#!/usr/bash

ID=07 # 01

for n in $(seq 0 31)  # set first and final model number 
    do
        echo 'Computing IGM attenuation for Lyn = ', ${n}
	
	python IGM_attenuation_lyman_series.py output/f_v12_${ID}.out output/IGM_attenuation_${ID}_LAF_ly${n}.output LAF ${n} >> log/IGM_attenuation_${ID}_LAF_ly${n}.log &
#	python IGM_attenuation_lyman_series.py output/f_v12_${ID}.out output/IGM_attenuation_${ID}_LLS_ly${n}.output LLS ${n} >> log/IGM_attenuation_${ID}_LLS_ly${n}.log &
#	python IGM_attenuation_lyman_series.py output/f_v12_${ID}.out output/IGM_attenuation_${ID}_DLA_ly${n}.output DLA ${n} >> log/IGM_attenuation_${ID}_DLA_ly${n}.log &

        echo 'waiting...'
	wait
    done
