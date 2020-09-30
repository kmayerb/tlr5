"""
Process Batch Files in resources directory
"""
import os 
resources = 'reads' #'/Volumes/Samsung_T5/kublin/106_fastq'
ft = [(f,os.stat(os.path.join(resources, f)).st_size) for f in os.listdir(resources) if f.endswith("fastq.gz")]
ft = sorted(ft,key= lambda x:x[1])
files_to_process = [x[0] for x in ft]

import time
for f in files_to_process:
	if os.path.isfile(f"outputs/{f}.flagellin.summary.tsv"):
		print(f"ALREADY PROCESSED : {f}")
	else:
		cmd = f"python bin/flagellin.py --input_fastq_gz {os.path.join(resources, f)} --outfile_filename outputs/{f}.flagellin.summary.tsv --threads 6"
		print(cmd)
		time.sleep(10)
		os.system(cmd)
		time.sleep(30)
