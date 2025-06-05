

This repository orgnizes new data release upon distirbution by merging Raw and Clean data
in the same sample directory and distribuitng it in their corresponding campaigns

Moreover, provides an easy to use script for uplaoding Metagnomics samples to our cluster for computational jobs

Create a directory and Download data in a location, Example:
```bash

mkdir data_release
wget -i X204SC24072989-Z02-F007.csv
unzip result_X204SC24072989-Z02-F007

```



Merge raw and clean data together as 02.Combined and distribute the samples in their respective location using campaign ID
- e Select either 16_S or 18_S or ITS based on the sample type data
- d newly downloaded file location
- f Parent final dir destination (NOTE: Directories are automatically selected based on the following google sheet )

```bash

python scripts/data_handling.py 
-e 16_S 
-d /data_release/X204SC24072989-Z02-F007/result_X204SC24072989-Z02-F007 
-f ../../sequencing_data/

```
GEOMOSAIC UPLOAD

The following script uploads forward and reverse files for each sample declared in a .txt file, along with a sample_table_{campaign}.tsv file used by GEOMOSAIC pipeline

```bash

python create_sample_table.py  
-d /media/edotacca/Loki/sequencing_data/Campaign/Metagenomes/ 
-s /sample_list.txt
-u 'user_name'

```
