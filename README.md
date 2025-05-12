This repository orgnizes new data release upon distirbution by merging Raw and Clean data
in the same sample directory and distribuitng it in their corresponding campaigns

STEP0
Create a directory and Download data in a location
```bash

mkdir data_release
wget -i data_release.csv
unzip result_*

```

STEP1

Merge raw and clean data together as 02.Combined and distribute the samples in their respective location using campaign ID

```bash

python scripts/data_distribution.py -e 18_S -d /data_release_21_02_Novo_gene/result_X204SC24072989-Z02-F007/ -f ../../sequencing_data/

```
The following script is used to merge raw data files for Geomosaic pipeline,
along for generating smaple_table.tsv

```bash

python create_sample_table.py  
-d /media/edotacca/Loki/sequencing_data/Campaign/Metagenomes/ 
-f /media/edotacca/Thor/geomosaic_upload/ 
-p 'G*'

```
