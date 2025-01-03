#!/bin/bash

#####Set Scheduler Configuration Directives#####
#Name the job:
#$ -N skm_libnorm_lbt005_selected_kwgrid_102024
#Send zhihan-chen@uiowa.edu at beginning/end/suspension of job
#$ -m bes

#E-mail address to send to
#$ -M zhihan-chen@uiowa.edu


#####End Set Scheduler Configuration Directives#####


#####Resource Selection Directives#####
#Select the queue to run in


#Request four cores
#$ -pe smp 36

#####End Resource Selection Directives#####

##### Analysis Syntax #####




module load stack/2021.1
module load py-tensorflow
module load stack/legacy
module load python/3.7.0

cd /old_Users/zchen190/soft\ kmeans/run\ file/102024
python run_argon.py
