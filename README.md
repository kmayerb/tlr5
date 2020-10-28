# tlr5

## Instructions

0. Ensure environment has a version of diamond and biopython

```bash
conda install -c bioconda diamond
pip install biopython
```

1. Clone this repository onto a working directory on the cluster. You will need the files present in the `inputs/` directory. 

2. A commandline program is available in `bin/flagellin.py`

* The primary input is a compressed fastq file of metagenomic reads.
* The output is a tabular summary of the best hits against a database of flagellin genes

```bash
python bin/flagellin.py --input_fastq_gz {os.path.join(resources, f)} --outfile_filename outputs/{f}.flagellin.summary.tsv --threads 6"
```

```bash
python bin/flagellin.py --help
usage: flagellin.py [-h] --input_fastq_gz INPUT_FASTQ_GZ --outfile_filename OUTFILE_FILENAME [--base_dir BASE_DIR] [--threads THREADS] [--tlr5_pattern TLR5_PATTERN] [--pfam_ncbi_fasta_filename PFAM_NCBI_FASTA_FILENAME] [--hu_accession_list_filename HU_ACCESSION_LIST_FILENAME] [--hu_table_filename HU_TABLE_FILENAME]
                    [--get_detail]

optional arguments:
  -h, --help            show this help message and exit
  --input_fastq_gz INPUT_FASTQ_GZ
                        Filename for an input fastq.gz file.
  --outfile_filename OUTFILE_FILENAME
                        Filename for output table summarizing flagellin discovered in the fastq.gz.
  --base_dir BASE_DIR   The base directory (OPTIONAL AND NOT NECESSARY)
  --threads THREADS     Max number of threads
  --tlr5_pattern TLR5_PATTERN
                        The tlr5_pattern examined in region 80-100 of the D1 sububit.
  --pfam_ncbi_fasta_filename PFAM_NCBI_FASTA_FILENAME
                        Filename for the PFAM NCBI protein fasts file (sequences should not have gaps).
  --hu_accession_list_filename HU_ACCESSION_LIST_FILENAME
                        Filename wiht list of NCBI accessions from Hu and Reeves 2020.
  --hu_table_filename HU_TABLE_FILENAME
                        Filename Table from Hu and Reeves 2020
  --get_detail          If this flag is used the full diamond tabular result is written in addition to the summary information
```


Calling this program will generate the following temporary files:
* flagellin_db.fasta
* flagellin_db.dmnd
* diamond.tabular_results.txt

However, the primary output will be found in outputs (e.g., outputs/129650059_TAAGGCGA-CTAAGCCT.fastq.gz.flagellin.summary.tsv). 

Each row corresponds with a single read

* `qseqid` - read ID, e.g, DH1DQQN1:1210:HKHNFBCX2:1:1101:11106:31171
* `sseqid` - best hit ID. e.g., ADF57343 
* `pident` - percent identity
* `bitscore` - bitsocre
* `Group` - Taxonomic Group
* `SuperKingdom` - Taxonomic SuperKingdom
* `Phylum` - Taxonomic Phylum
* `Class` - Taxonomic Class
* `Order` - Taxonomic Order
* `Family` - Taxonomic Family
* `Genus` - Taxonomic Genus
* `Species`- Taxonomic Species
* `pfid` - descriptive id of best hist (e.g., Roseburia inulinivorans	ADF57343/3-139)
* `seq` - query sequencce
* `re_search` -  e.g., LQRMNEL
* `tlr5_status` - TLR5 activating status of all hits above threshold [TLR5+ , TRL5-, or TLR5+/-] 
* `percent_positive` - [0-1] percentage of hits that are TRL5+
* `n` - number of hits


3. The commandline script may be called on many fastq files within a target folder as shown 
in `run.py`
