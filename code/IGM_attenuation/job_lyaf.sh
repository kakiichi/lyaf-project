#!/usr/bin/tcsh

#$ -e run.err
#$ -o run.out

#$ -l h_vmem=1000M    
#$ -cwd
#$ -m be
#$ -M kakiichi@mpa-garching.mpg.de
#$ -l hostname="lnx-1[12397]"

unlimit 

# DLA regime
python tau_eff_v1.0.py output/f_v12_${ID}.out output/tau_eff_DLA_${ID}.output DLA >> log/tau_eff_DLA_${ID}.log 

# LLS regime
#python tau_eff_v1.0.py output/f_v12_${ID}.out output/tau_eff_LLS_${ID}.output LLS >> log/tau_eff_LLS_${ID}.log 

# LAF regime
#python tau_eff_v1.0.py output/f_v12_${ID}.out output/tau_eff_LAF_${ID}.output LAF >> log/tau_eff_LAF_${ID}.log 
