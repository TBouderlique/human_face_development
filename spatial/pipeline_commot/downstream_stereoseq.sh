#!/bin/bash
batch=C01939B2
db="Secreted Signaling"


for pathway in PATHWAY1 PATHWAY2 PATHWAY3
do
   echo $pathway
   source /home/felix/miniconda3/bin/activate commot
   taskset -c 80,81,82 python downstream_deg.py $batch $pathway $db

done
echo end