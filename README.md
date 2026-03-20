

This repository cotains three scripts for easing common practices in day-t-day routine. Uploading of raw readds to a server to be ready for GEOMOSAIC setup, collection and gathering of basic statistics regarding results obtained after running geomosaic and distribution of newly released sequences from a common folder based on our sample_database.


### CREATE CONDA ENV
```bash
conda env create -f environment.yml
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

scp scripts/geomosaic_statistics.py {username}@ibiscohpc-ui.scope.unina.it:/{remote_path}

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
Either, pass a sample_table.tsv file such as:

| r1                         | r2                          | sample       |
| ---------------------------|-----------------------------|--------------|
| AC_280625_F_RBKMV.fastq.gz | AC_280625_F_RBKTMV.fastq.gz | AC_280625_F  |
| BC_200625_S_RBKMV.fastq.gz | BC_200625_S_RBKTMV.fastq.gz | BC_200625_S  |
| LS_230625_F_RBKMV.fastq.gz | LS_230625_F_RBKTMV.fastq.gz | LS_230625_F  |
| SF_221019_F_RBKMV.fastq.gz | SF_221019_F_RBKTMV.fastq.gz | SF_221019_F  |
| SF_221019_S_RBKMV.fastq.gz | SF_221019_S_RBKTMV.fastq.gz | SF_221019_S  |


Run the following command, substituting with the right paths and flags:
```bash

python scripts/geomosaic_uploader.py -h

usage: Upload raw seqs AND sample_table.tsv to a remote server [-h] [-d DATA_DIR] [-t SAMPLE_TABLE] [-c CAMPAIGN_NAME] [-p PATTERN] [-s SUBSET_SAMPLES] [-u IBISCO_USER] [-z]

options:
  -h, --help            show this help message and exit
  -d DATA_DIR, --data_dir DATA_DIR
                        Provide directory ABSOLUTE path storing sample sequences
  -t SAMPLE_TABLE, --sample_table SAMPLE_TABLE
                        Provide geomosaic ABSOLUTE path to sample table
  -H REMOTE_HOST, --remote_host REMOTE_HOST
                        Provide remote server hostname or IP (e.g. dgiovannelli@ibiscohpc-ui.scope.unina.it)
  -p REMOTE_PATH, --remote_path REMOTE_PATH
                        Provide absolute remote path on the server (e.g. /ibiscostorage/GiovannelliLab/raw)
  -c CAMPAIGN_NAME, --campaign_name CAMPAIGN_NAME
                        Provide campaign name: Ex ARG23 or TEST
  -z, --dry_run         Will attempt a dry run without uploading files
```
Example:
```bash

```

### DATA-HANDLING

DOWNLOAD DATA RELEASE
Create a directory and Download data in a location, Example:
```bash

mkdir data_release
wget -i X204SC24072989-Z02-F007.csv
unzip result_X204SC24072989-Z02-F007

```
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