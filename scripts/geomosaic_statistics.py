import os
import glob
import shutil
import argparse
from Bio import SeqIO
import pandas as pd
import re
from tqdm import tqdm
from contextlib import redirect_stdout
import json
import sys

####NOTA####
# this script is executed by typing on terminal:
# conda activate geomosaic
# python geomosaic_statistics.py -w $campaign/geomosaic/ -s $campaign/list_samples.txt 


def main():

    args = parse_args()

    samples, output_dir = get_samples(
            geo_samples=args.samples,
            output_dir=args.output_dir
    )

    compute_statistics(
            geomosaic_work_dir=args.working_dir,
            output_dir=output_dir,
            sample_list=samples
    )


def get_samples(geo_samples:str,output_dir:str):
    
    os.makedirs(output_dir, exist_ok=True)

    if os.path.exists(geo_samples):
        
        try:
            with open(geo_samples,'r') as reader:
                samples = [line.strip() for line in reader]
                if samples:
                    return sorted(samples), output_dir
                else:
                    print("Sample list is empty")
                    sys.exit()
    
        except Exception as error:
            return None

    else:
        print(f'No files named: {geo_samples} found')
        return []

def check_exists(file:str):

    if os.path.exists(file):

        try:
            dataframe = pd.read_csv(file, sep=',')
            return dataframe
            
        except Exception as e:
            print(f"Error while reading {file}")
            return None


def get_tax_classified_reads(geomosaic_work_dir:str,samples:list,pckg:str,output:str):
    
    taxa_rank = "genus.tsv"
    taxa = taxa_rank.split('.')[0]
    file_output = os.path.join(output,f'kaiju_{taxa}_count.csv')

    dataframe = check_exists(file_output)

    if dataframe is not None:
        print('[+] Loading previous KAIJU report')
        return dataframe, file_output
    
    else:
        all_counts = []
        print('[+] Parsing KAIJU output tables')

        for s in tqdm(samples, desc="Processing taxa classification"):
            pckg_dir = os.path.join(geomosaic_work_dir,s,pckg)

            file = os.path.join(pckg_dir,taxa_rank)
            if os.path.exists(file):
                # this give paired-reads NOT single reads count
                dataframe = pd.read_csv(file,sep="\t")
                unclassified = ["unclassified",f"cannot be assigned to a (non-viral) {taxa}"]
                unclass_reads = dataframe.loc[dataframe["taxon_name"].isin(unclassified), "reads"].sum()

                taxa_table = dataframe.set_index("taxon_name")
                taxa_table.drop(labels = unclassified,axis = 0, inplace=True)
                class_reads = taxa_table["reads"].sum()

                all_counts.append({'sample': s,'class_read-pairs': int(class_reads),'unclass_read-pairs' : int(unclass_reads)})
            else:
                print(f"File {file} NOt found")
                continue
    
    dataframe = pd.DataFrame(all_counts, columns = ['sample','class_read-pairs','unclass_read-pairs'])
    dataframe.to_csv(file_output, sep = ',', index = False)

    return dataframe, file_output
                    


def get_fastp_stats(geomosaic_work_dir:str,samples:list,pckg:str,output:str):

    file_output = os.path.join(output,'fastp_count.csv')
    dataframe = check_exists(file_output)

    if dataframe is not None:
        print('[+] Loading previous FASTP report')
        return dataframe, file_output
    
    else:
        print('[+] Parsing FASTP report')
        all_samples = []
        for s in tqdm(samples, desc="Processing reads"):
            
            pckg_dir = os.path.join(geomosaic_work_dir,s,pckg)
            if os.path.exists(pckg_dir):
                file = os.path.join(pckg_dir,'report.json')

                if os.path.exists(file):
                    with open(file, "r") as f:
                        data = json.load(f)

                        summary = data["summary"]
                        filtering = data["filtering_result"]

                        stats = {
                        "sample": s,
                        "total_reads_raw": summary['before_filtering']['total_reads'],
                        "total_reads_clean": summary['after_filtering']['total_reads'],
                        "discarded_reads" : int(summary['before_filtering']['total_reads'] - summary['after_filtering']['total_reads']),
                        "gc_content": round(summary['after_filtering']['gc_content'] * 100,2),
                        "duplication": round(data['duplication']['rate'] * 100,2)
                        }

                        all_samples.append(stats)
            else:
                print(f'No {pckg_dir} found')
                continue
    
    dataframe = pd.DataFrame(all_samples, columns = ['sample','total_reads_raw','total_reads_clean','discarded_reads','gc_content','duplication'])
    dataframe.to_csv(file_output, sep = ',', index = False)
    
    return dataframe, file_output



def get_fastp_readcounts(geomosaic_work_dir:str,samples:list,pckg:str,output:str):

    file_output = os.path.join(output,'reads_count.csv')
    dataframe = check_exists(file_output)

    if dataframe is not None:
        print('[+] Loading previous READSQC report')

        return dataframe, file_output
    
    else:
        print('[+] Parsing READSQC report')
        all_counts = []
        
        for s in tqdm(samples, desc="Processing reads"):
            
            pckg_dir = os.path.join(geomosaic_work_dir,s,pckg)
            if os.path.exists(pckg_dir):

                file = os.path.join(pckg_dir,'geomosaic_readscount.txt')
                if os.path.exists(file):

                    with open(file,'r') as reader:
                        count = reader.readline().strip()
                        all_counts.append({'sample': s, 'reads': int(count)})
                else:
                    print(f'No {file} found')
            else:
                print(f'No {pckg_dir} found')
                continue

        dataframe = pd.DataFrame(all_counts, columns = ['sample','reads'])
        dataframe.to_csv(file_output, sep = ',', index = False)

        return dataframe, file_output


def get_contigs_length(geomosaic_work_dir:str,samples:list,pckg:str,output):

    #HEADER_PATTERN = re.compile(r'flag=(\d+)\s+multi=([\d\.]+)\s+len=(\d+)')
    file_output = os.path.join(output,'contigs_stats.csv')
    dataframe = check_exists(file_output)

    if dataframe is not None:
        return dataframe, file_output

    else:
        print('[+] Parsing assembly_contigs..')
        all_contigs_lengths= []

        for s in tqdm(samples, desc="Processing contigs"):
            pckg_dir = os.path.join(geomosaic_work_dir,s,pckg)

            total_assembly_length = 0
            n_contigs = 0

            if os.path.exists(pckg_dir):
                file = os.path.join(pckg_dir,'filtered_contigs.fasta')
                
                try:
                        for record in (SeqIO.parse(file,'fasta')):
                            total_assembly_length += len(record.seq)
                            n_contigs += 1
                        all_contigs_lengths.append({'sample': s, 'assembly_lenght': total_assembly_length, 'n_contigs': n_contigs})
                            
                except Exception as e:
                    print(f"An error occurred during file parsing: {e}")
                    return

        dataframe = pd.DataFrame(all_contigs_lengths, columns = ['sample', 'assembly_lenght','n_contigs'])
        dataframe.to_csv(os.path.join(file_output), sep = ',', index = False)

    return dataframe, file_output



def get_mags_stats(source_path:str,samples:list,pckg:str,output):

    all_dataframes = []
    
    main_file_out = os.path.join(output,'all_samples_mags.tsv')
    dataframe = check_exists(main_file_out)

    if dataframe is not None:
        return dataframe, main_file_out
    
    print('[+] parsing MAGs stats file..')
    for s in tqdm(samples, desc="Processing MAGs"):

        pckg_dir = os.path.join(source_path,s,pckg)
        file_out = os.path.join(output,f'{s}_MAGs.tsv')
        
        if os.path.exists(pckg_dir):
            file = os.path.join(pckg_dir,'MAGs.tsv')
            # if os.path.exists(file) and os.path.exists(file_out):
            #     return 
            try:
                df_checkm = pd.read_csv(file, sep='\t')
            except FileNotFoundError:
                tqdm.write(f"CheckM file not found for sample {s}. Skipping.")
                continue

            mags = df_checkm['MAGs'].to_list()
            all_mags_lengths = []

            for mag in mags:
                
                total_mag_lenght = 0
                n_contigs = 0

                MAG = os.path.join(pckg_dir,'fasta',f'{mag}.fa')
                if os.path.exists(MAG):
                    
                    for record in (SeqIO.parse(MAG,'fasta')):
                        total_mag_lenght += len(record.seq)
                        n_contigs += 1
                        
                mag_row =df_checkm[df_checkm['MAGs'] == mag]

                complet = mag_row['Completeness'].values[0]
                contam = mag_row['Contamination'].values[0]

                weighted_complet = (complet/100)*total_mag_lenght
                weighted_contam = (contam/100)*total_mag_lenght

                all_mags_lengths.append({'MAGs': mag, 'sample' : s, 'Completeness': complet, 'Contamination': contam, 
                                         'weighted_completness': weighted_complet, 'weighted_contamination': weighted_contam,
                                          'total_length_bp': total_mag_lenght, 'n_contigs': n_contigs})

            df = pd.DataFrame(all_mags_lengths,columns = ['MAGs','sample','Completeness','Contamination',
                                                    'weighted_completeness','weighted_contamination',
                                                    'total_length_bp','n_contigs'])
        
        df.to_csv(file_out, sep = '\t', index = False)
        all_dataframes.append(df)                                                  

    dataframe = pd.concat(all_dataframes, ignore_index=True)
    dataframe.to_csv(main_file_out, sep = '\t', index = False)
    
    return dataframe, main_file_out



def average(dataframe:pd.DataFrame,var:str):
    
    if var in dataframe.columns and not dataframe[var].empty:

        if var in ['Completeness','Contamination']:
            avg = dataframe[var].mean()
            return round(avg,2)
        else:
            avg = dataframe[var].mean()
            return round(avg) 
    else:
        return 0
        

def compute_statistics(geomosaic_work_dir:str, output_dir:str, sample_list:list,output_log_file='shell_stats.txt'):

    if sample_list:

        df_tax_stats, file_kaiju = get_tax_classified_reads(geomosaic_work_dir,sample_list,'kaiju', output_dir)

        df_fastp_stats, file_fastp = get_fastp_stats(geomosaic_work_dir, sample_list, 'fastp', output_dir)

        df_reads_counts,file_reads = get_fastp_readcounts( geomosaic_work_dir, sample_list, 'fastqc_readscount', output_dir)

        df_assembly_lengths,file_contigs = get_contigs_length( geomosaic_work_dir, sample_list, 'megahit', output_dir)

        #df_mags_stats , file_mags = get_mags_stats(geomosaic_work_dir,sample_list,'mags',output_dir)
        print('\n')

        all_stats_df = [df_fastp_stats, df_reads_counts, df_tax_stats, df_assembly_lengths]

        # gather reads and assembly based stats
        
        stats_table = pd.DataFrame(sorted(sample_list), columns=["sample"])
        
        for df in all_stats_df:
            if df is not None:
                temp = pd.merge(stats_table,df ,how='left', on='sample')
                stats_table = temp.copy()
        
        print('\n',stats_table,'\n')
        table_output = os.path.join(output_dir,"stats_table.tsv")
        stats_table.to_csv(table_output, sep="\t", index = False)

        '''
        
        log_file = os.path.join(output_dir,output_log_file)
        with open(log_file, 'wt') as f:
            
            with redirect_stdout(f):

        #----------------------------------------------------------------------------------------------------
                # calcualte average reads per sample:
                if df_reads_counts.empty:
                    df_reads_counts = pd.read_csv(file_reads,sep = ',')

                min_len = df_reads_counts['reads'].min()
                max_len = df_reads_counts['reads'].max()
                print(f"Total reads length range (bp): {min_len:,} to {max_len:,}")

                avg_total_reads = average(df_reads_counts, 'reads')
                print(f"Average total read length (bp): {avg_total_reads/1_000_000:.2f} million \n")
                #----------------------------------------------------------------------------------------------------
                # Calculate contigs range
                if df_assembly_lengths.empty:
                    df_assembly_lengths = pd.read_csv(file_contigs,sep = ',')
                    
                min_len = df_assembly_lengths['total_length_bp'].min()
                max_len = df_assembly_lengths['total_length_bp'].max()
                print(f"Total assembly length range (bp): {min_len:,} to {max_len:,}")

                avg_total_length = average(df_assembly_lengths, 'total_length_bp')
                print(f"Average total assembly length (bp): {avg_total_length/1_000_000:.2f} million \n")        

                #----------------------------------------------------------------------------------------------------
                #this function gather stats result from reads,contigs and mags in na unique (per samaple averaged) way
                if df_mags_stats.empty:
                    df_mags_stats = pd.read_csv(file_mags,sep = '\t')

                avg_total_completness = average(df_mags_stats, 'Completeness')
                avg_total_contamination = average(df_mags_stats, 'Contamination')

                print(f"Average total completness (unweighted): {avg_total_completness} \n")
                print(f"Average total contamination (unweighted): {avg_total_contamination} \n")
            
            print(f"\n[+] Statistics log written to: {output_log_file}")
    else:
        return None

'''

def parse_args():
    parser = argparse.ArgumentParser("TEXT HERE")
    
    parser.add_argument("-w", "--working_dir",
        help="Absolute apth to geomosaiuc working directory.",
        type=str,
        required=True
    )
    parser.add_argument("-s", "--samples",
        help="Samples file list.",
        type=str,
        required=True
    )
    parser.add_argument("-o", "--output_dir",
        help="Where the files are wirtten.",
        type=str,
        required=True
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()