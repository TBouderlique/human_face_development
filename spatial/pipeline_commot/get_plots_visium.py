import sys
import os
samp = sys.argv[1]
os.system(f"python get_pathway_plots_visium.py {samp}")
os.system(f"python get_sender_reciever_plots_visium.py {samp}")
