import sys
import os
samp = sys.argv[1]


os.system(f"python get_pathway_plots_stereoseq.py {samp}")
os.system(f"python get_sender_reciever_plots_stereoseq.py {samp}")
