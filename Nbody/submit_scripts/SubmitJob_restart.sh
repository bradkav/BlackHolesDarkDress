#!/bin/bash

cd /home/kavanagh/PBH/

./PrepareJob_restart.sh $1 $2 $3

#echo "Submitting Job: RunPBH_$1.sh"
cd /home/kavanagh/PBH/sims/$1
qsub RunPBH_$1.sh
