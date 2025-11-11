#!/bin/bash
# 070623_A1_eye_socket_visium
# 070623_B1_sagital-mandible
# 070623_D1_W8_frontal_mandible_visium


batch=070623_D1_W8_frontal_mandible_visium
db=secreted_signaling
for pathway in ncWNT PTN MK ANGPT CALCR VISFATIN ANGPTL MIF CXCL GDF FGF PERIOSTIN PDGF TGFb GAS SEMA3 PROS SPP1 WNT BMP VEGF IGF
do
   echo $pathway
   source /home/felix/miniconda3/bin/activate commot
   taskset -c 2,3 python downstream_deg.py $batch $pathway $db
done
echo end