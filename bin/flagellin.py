#!/usr/bin/env python3

"""
September 2020
kmayerbl@fredhutch.org

flagellin.py

formerly: pipeline_flagellin_TLR5.py

The D1 subunit of bacterial flagellin has the potential to act an an immune adjuvant via 
direct iteraction with Toll LIke Receptor 5. Structural studies suggest an Arginine at reference position 
R[89] in bascillus subtilis is conserved amoung TLR5-activiating flagelliln. 

Based on multiple sequence alignment, the molecular pattern used in this study t
o identify TLR5 activating (hereafter TLR5+) flagelin is '([LVIM].R[MIALV]..[LI])'.

This pipeline permits the interogation of raw metagenomic .fastq datafiles, quantifying detectable bacterial flagellins
and scoring whether they match TLR5+ reference genes. 

The primary motif responsible for the D1-flagelline subunit-TLR5 interface is found within 
the Bacterial flagellin N-terminal helical region (https://pfam.xfam.org/family/PF00669). 

As of September 2020, EMBL-EBI PFAM database conatins (73m795) protein sequences from NCBI annotated
as members of thsi protein family. These secquencs were downloaded /data/PF00669_ncbi.txt/.

He et al 2020 also constructed a taxanomically annotated database of flagellin sequences from the 
Protein Family00669 as well Protein Family 00700 Bacterial flagellin C-terminal helical region).

Version 1 of our flagellin database is the intersection of sequences in PF00669_ncbi.txt and 
the Flagellin DB. 
"""

import numpy as np
import os
import pandas as pd 
import re

from  Bio import SeqIO

def safe_re_search(s, pattern, group = 0):
	r = re.search(pattern =pattern , string =s)
	if r is not None:
		return r.groups()[group]
	else:
		return None

def construct_reference_flagellin_dataframe(pfam_ncbi_fasta_filename,
		hu_accession_list_filename,
		hu_table_filename,
		tlr5_pos_pattern = '([LVIM].R[MIALV]..[LI])',
		base_dir = '/Users/kmayerbl/TLR/'):

	"""
	This dataframe will be used later on to 
	score diamond results
	"""
	record_dict = SeqIO.to_dict(SeqIO.parse(pfam_ncbi_fasta_filename, 'fasta'))

	mapr = {x.strip().split("/")[0]:x for x in record_dict.keys()}
	ncbi_in_pfam = [x.strip().split("/")[0] for x in record_dict.keys()]

	with open(hu_accession_list_filename, 'r') as f:
		accessions = f.readlines()
	
	hu_ncbi_all = [x.strip() for x in accessions]
	set(hu_ncbi_all).intersection(set(ncbi_in_pfam))
	acc = list(set(hu_ncbi_all).intersection(set(ncbi_in_pfam)))
	df = pd.DataFrame({"Accession": acc})


	hu = pd.read_csv(hu_table_filename)
	df = df.merge(hu, how ="left") 
	df['pfid']   = df["Accession"].apply(lambda x: mapr.get(x))
	df['seq']    = df['pfid'].apply(lambda x: str(record_dict[x].seq))
	df['seq_region']   = df['seq'].apply(lambda x: x[80:100])
	
	df['re_search'] = df['seq_region'].apply(lambda x : safe_re_search(x, pattern = tlr5_pos_pattern) )

	return df

def write_flaggelin_db_fasta(df, fasta_name = "flagellin_db.fasta"):
	with open(fasta_name, 'w') as oh:
		for i,r in df.iterrows():
			oh.write(f">{r['Accession']}\n{r['seq']}\n")
	return fasta_name


def make_diamond_db(fasta_name ='flagellin_db.fasta', db_name = 'flagellin_db.dmnd'):
	cmd = f"diamond makedb --in {fasta_name} --db {db_name} --threads 2"
	print(cmd)
	os.system(cmd)
	return db_name 

def diamond_blastx(db_name = 'flagellin_db.dmnd', 
				   query_fastq_gz = 'Volumes/Samsung_T5/kmayerbl/cf_input/C94BLACXX_4_GCTACGCT_GTAAGGAG.S.fq.fastq.gz',
				   diamond_output_filename = 'diamond.tabular_results.txt',
				   min_percent_identity = 75,
				   output_format = 6,
				   threads = 6):
	cmd = f"diamond blastx --db {db_name} --query {query_fastq_gz }	--out {diamond_output_filename} --outfmt {output_format} --id {min_percent_identity} --threads {threads}"   
	print(cmd)
	os.system(cmd)
	return diamond_output_filename

def adjudicate_list(l):
	scores = [1 if (x is not None) else 0 for x in l]
	percent_positive = np.sum(scores) / len(scores)
	
	if percent_positive == 1.0:
		tlr5_status = "TLR5+"
	
	elif percent_positive > 0.0:
		tlr5_status  = "TRL5+/-"
	
	elif percent_positive == 0.0:
		tlr5_status  = "TLR5-"
	
	return tlr5_status , percent_positive, len(scores)


def parse_diamond_blastx_output(df, diamond_output_filename = 'diamond.tabular_results.txt'):

	diamond_df = pd.read_csv( diamond_output_filename, sep = "\t", header = None)
	diamond_df .columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
	diamond_df  = diamond_df.merge(df, how ="left", left_on = "sseqid", right_on = "Accession")

	summary_lines = list()

	for i,g in diamond_df.groupby(['qseqid']):

		max_bitscore = g['bitscore'].max()
		gmax = g[g['bitscore'] == max_bitscore].copy()
		gmax[['qseqid','sseqid','pident','bitscore','Group', 'SuperKingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species', 'pfid', 'seq', 're_search']]
		TLR5_interfaces = gmax['re_search'].to_list()
		tlr5_status , percent_positive, n  = adjudicate_list(TLR5_interfaces)
		first_row = gmax[['qseqid','sseqid','pident','bitscore','Group', 'SuperKingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species', 'pfid', 'seq', 're_search']].iloc[0,:].copy()
		first_row['tlr5_status']      = tlr5_status
		first_row['percent_positive'] = percent_positive
		first_row['n']                = n
		summary_lines.append(first_row)        
	
	summary_df = pd.DataFrame(summary_lines)      
	
	return summary_df, diamond_df 


if __name__ == "__main__":
	""" 
	python bin/flagellin.py \
	--input_fastq_gz /Volumes/Samsung_T5/kmayerbl/cf_input/C94BLACXX_4_GCTACGCT_GTAAGGAG.S.fq.fastq.gz  \
	--outfile_filename outputs/flagellin_outputfile.tsv \
	--get_detail
	"""
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('--input_fastq_gz',
		type = str, action="store", required=True, help ="Filename for an input fastq.gz file.")
	parser.add_argument('--outfile_filename', 
		action="store", type = str, required=True, help = "Filename for output table summarizing flagellin discovered in the fastq.gz.")
	parser.add_argument('--base_dir',
		type = str, action="store", required=False, default = '', help ="The base directory (OPTIONAL AND NOT NECESSARY)")
	parser.add_argument('--tlr5_pattern',
		type = str, action="store", required=False, default = '([LVIM].R[MIALV]..[LI])', help ="The tlr5_pattern examined in region 80-100 of the D1 sububit.")
	parser.add_argument('--pfam_ncbi_fasta_filename', 
		action="store", type = str, required=False, default = 'inputs/PF00669_ncbi.txt', help = "Filename for the PFAM NCBI protein fasts file (sequences should not have gaps).")
	parser.add_argument('--hu_accession_list_filename', 
		action="store", type = str, required=False, default = 'inputs/hu_all_ncbi.txt', help = "Filename wiht list of NCBI accessions from Hu and Reeves 2020.")
	parser.add_argument('--hu_table_filename', 
		action="store", type = str, required=False, default = 'inputs/hu_table.csv', help = "Filename Table from Hu and Reeves 2020")
	parser.add_argument('--get_detail', 
		action='store_true', default = False, help = "If this flag is used the full diamond tabular result is written in addition to the summary information")
	

	args = parser.parse_args()

	for arg in vars(args):
	    print(f"{arg}={getattr(args, arg)}")
	
	input_fastq_gz              = args.input_fastq_gz 

	input_fastq_gz              = args.input_fastq_gz	            #'/Volumes/Samsung_T5/kmayerbl/cf_input/C94BLACXX_4_GCTACGCT_GTAAGGAG.S.fq.fastq.gz' 
	tlr5_pattern                = args.tlr5_pattern  				# '([LVIM].R[MIALV]..[LI])'
	base_dir                    = args.base_dir 					# '/Users/kmayerbl/TLR/'
	pfam_ncbi_fasta_filename    = args.pfam_ncbi_fasta_filename  	# os.path.join(base_dir, 'inputs','PF00669_ncbi.txt')
	hu_accession_list_filename  = args.hu_accession_list_filename 	# os.path.join(base_dir,'inputs', 'hu_all_ncbi.txt')
	hu_table_filename 			= args.hu_table_filename			# os.path.join(base_dir, 'inputs', 'hu_table.csv')
	outfile_filename            = args.outfile_filename             # "test.outfile.txt"

	ref_df = construct_reference_flagellin_dataframe(pfam_ncbi_fasta_filename= pfam_ncbi_fasta_filename,
		hu_accession_list_filename = hu_accession_list_filename, 
		hu_table_filename = hu_table_filename,
		tlr5_pos_pattern =tlr5_pattern, 
		base_dir = '/Users/kmayerbl/TLR/')
	
	print(f"REFERENCE DATAFRAME:")
	print(ref_df)

	fasta_name      = write_flaggelin_db_fasta(df =ref_df, fasta_name = "flagellin_db.fasta")
	print(f"WROTE {fasta_name}")
	
	db_name         = make_diamond_db(fasta_name ='flagellin_db.fasta', db_name = 'flagellin_db.dmnd')
	print(f"WROTE {db_name}")
	
	output_filename = diamond_blastx(  db_name = db_name, 
									   query_fastq_gz = input_fastq_gz,
									   diamond_output_filename = 'diamond.tabular_results.txt',
									   min_percent_identity = 75,	
									   output_format = 6,
									   threads = 6)
	print(f"WROTE {output_filename}")

	summary_df, diamond_df = parse_diamond_blastx_output(diamond_output_filename = 'diamond.tabular_results.txt', df = ref_df)

	summary_df.to_csv( outfile_filename, sep = "\t") 
	if args.get_detail:
		print(f"{outfile_filename}.full_diamond.blastx.tsv")
		diamond_df.to_csv( f"{outfile_filename}.full_diamond.blastx.tsv", "\t")