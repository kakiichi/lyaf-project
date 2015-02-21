#!/usr/bin/tcsh 


for ID in {01..10}  # set first and final model number 
    do
        echo 'Submit the job ID ...',${ID}      

        qsub -v ID=$ID job_lyaf.sh

        echo ''
    done
