

This repository orgnizes new data release upon distirbution by merging Raw and Clean data
in the same sample directory and distribuitng it in their corresponding campaigns
Moreover, provides an easy to use script for uplaoding Metagnomics samples to our cluster for computational jobs

DOWNLOAD DATA RELEASE
Create a directory and Download data in a location, Example:
```bash

mkdir data_release
wget -i X204SC24072989-Z02-F007.csv
unzip result_X204SC24072989-Z02-F007

```

### CREATE CONDA ENV
```bash
conda env create -f environment.yml
```
### MERGE AND DISTRIBUTE FILES

Merge raw and clean data together as 02.Combined and distribute the samples in their respective location using campaign ID
- -e Select either 16_S or 18_S or ITS based on the sample type data
- -d newly downloaded file location
- -f Parent final dir destination (NOTE: Directories are automatically selected based on the following google sheet )

```bash

python scripts/data_handling.py 
-e 16_S 
-d /data_release/X204SC24072989-Z02-F007/result_X204SC24072989-Z02-F007 
-f ../../sequencing_data/

```
### GEOMOSAIC STATISTICS

This script collects and gather statistics obtained running the GEOMOSAIC pipeline, in particular it gets you:
- fastp counts (raw and cleaned)
- contigs counts
- assembly stats -> assembly_stats.tsv
- mags statistics
For getting the information you want, you must ensure to have run the following modules:
- pre_processing, assembly, assembly_qc, binning_qa, and mags_retrieval modules
This script MUST BE uploaded to the server in your own directory, for example in your own user dir 
ex. user/etaccaliti

Upload the script to remote server:
```bash

scp scripts/geomosaic_statistics.py {username}@ibisco-hpc.ui.unina.it:/{remote_path}

```
Run the following command, substituting with the right paths and flags

```bash

python scripts/geomosaic_statistics.py -h
usage: Collect statistics and simple average [-h] -w WORKING_DIR -s SAMPLES -o OUTPUT_DIR

options:
  -h, --help            show this help message and exit
  -w WORKING_DIR, --working_dir WORKING_DIR
                        Absolute apth to geomosaiuc working directory.
  -s SAMPLES, --samples SAMPLES
                        Samples file list.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Where the files are wirtten.

```

### GEOMOSAIC UPLOADER

The following script uploads forward and reverse files for each sample declared in a sample_table_{campaign}.tsv file (The likes of geomosaic ). 

Run the following command, substituting with the right paths and flags:
```bash

python scripts/geomosaic_uploader.py -h
usage: Upload raw seqs AND sample_table.tsv to IBISCO server [-h] [-d DATA_DIR] [-t SAMPLE_TABLE] [-c CAMPAIGN_NAME] [-p PATTERN] [-s SUBSET_SAMPLES] [-u IBISCO_USER] [-z]

options:
  -h, --help            show this help message and exit
  -d, --data_dir DATA_DIR
                        directory storing sample sequences
  -t, --sample_table SAMPLE_TABLE
                        provides geomosaic sample table
  -c, --campaign_name CAMPAIGN_NAME
                        Provide campaign name: Ex ARG23 or TEST
  -p, --pattern PATTERN
                        Sequence pattern to search, default: *fq.gz
  -s, --subset_samples SUBSET_SAMPLES
                        provide a text file with each line a sample
  -u, --ibisco_user IBISCO_USER
                        Provide user account
  -z, --dry_run         Will attempt a dry run without uploading files

```
