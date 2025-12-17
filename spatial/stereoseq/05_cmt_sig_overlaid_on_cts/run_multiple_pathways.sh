#!/bin/bash
# Select pathways and which id of sample to run on

source /home/felix/miniconda3/bin/activate commot

for pathway in ncWNT PDGF MK PERIOSTIN BMP NRG IGF SEMA3 PTN
do
   echo $pathway
   python plt_commot_signals_overlaid_on_celltypes.py $pathway 7

done
echo end