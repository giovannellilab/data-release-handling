### AUTHOR: EDOARDO TACCALITI ###
## Date: April 2026 ###
## edoardotaccaliti@gmail.com ###

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
# import from my previous script 
from geomosaic_statistics import average

# This script is used to gather the MAGs from the different campaigns and put them in a single folder, with a single metadata file.
# The script also checks for the presence of the MAGs and their completeness and contamination, and filters them based on the provided thresholds.
def get_samples(geo_samples:str)->list:
    
    if os.path.exists(geo_samples):

        with open(geo_samples,'r') as reader:
            samples = [line.strip() for line in reader]
            return samples
    else:
        print(f'No files named: {geo_samples} found')
        return []


def copy_rename_mag(mag_file, output_dir, sample_id, glab_signature):
    
    mag_name = os.path.basename(mag_file).split('.')[0]  # Get the base name without extension
    extension = os.path.basename(mag_file).split('.')[1]  # Get the extension

    new_mag_name = f"{sample_id}_{glab_signature}_{mag_name}.{extension}"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(os.path.join(output_dir, new_mag_name)):
        shutil.copy(mag_file, os.path.join(output_dir, new_mag_name))
    else:
        tqdm.write(f"File {new_mag_name} already exists in {output_dir}. Skipping copy.")
    
    return new_mag_name


def collect_mags(s, pckg, working_dir, output_dir):
    """
    Collect MAG data for a single sample.
    Returns a DataFrame or None if the sample should be skipped.
    """

    pckg_dir = os.path.join(working_dir, s, pckg)
    file_out = os.path.join(output_dir, f'{s}_MAGs.tsv')

    if not os.path.exists(pckg_dir):
        tqdm.write(f'No {s} found at {pckg_dir}. Skipping.')
        return None

    file = os.path.join(pckg_dir, 'MAGs.tsv')
    
    try:
        df_checkm = pd.read_csv(file, sep='\t')
    except FileNotFoundError:
        tqdm.write(f"Error opening MAGs.tsv: file not found for sample {s}. Skipping.")
        return None

    mags = df_checkm['MAGs'].to_list()
    all_mags_lengths = []

    for mag in mags:

        total_mag_length, n_contigs = 0, 0

        MAG = os.path.join(pckg_dir, 'fasta', f'{mag}.fa')
        if not os.path.exists(MAG):
            tqdm.write(f"FASTA file for MAG {mag} not found in sample {s}. Skipping.")
            continue

        for record in SeqIO.parse(MAG, 'fasta'):
            total_mag_length += len(record.seq)
            n_contigs += 1

        #renaming the MAG file and copying it to the output directory with a new name that includes the sample id and a signature for GLabCoEvo
        glab_mag = copy_rename_mag(MAG, output_dir=output_dir, sample_id=s, glab_signature='GlabCoEvo')

        mag_row = df_checkm[df_checkm['MAGs'] == mag]
        complet = mag_row['Completeness'].values[0]
        contam = mag_row['Contamination'].values[0]

        weighted_complet = (complet / 100) * total_mag_length
        weighted_contam = (contam / 100) * total_mag_length

        row = {
            'GLab_mag': glab_mag,
            'MAGs': mag,
            'binID': mag_row['binID'].values[0],
            'sample': s,

            'Marker lineage': mag_row['Marker lineage'].values[0],
            'n_genomes': mag_row['# genomes'].values[0],
            'n_markers': mag_row['# markers'].values[0],
            'n_marker_sets': mag_row['# marker sets'].values[0],
            'Completeness': complet,
            'Contamination': contam,
            'Strain heterogeneity': mag_row['Strain heterogeneity'].values[0],

            'weighted_completeness': weighted_complet,
            'weighted_contamination': weighted_contam,
            'total_length_bp': total_mag_length,
            'n_contigs': n_contigs
        }
        all_mags_lengths.append(row)

    df_samples_mags = pd.DataFrame(all_mags_lengths)
    #df_samples_mags.to_csv(file_out, sep='\t', index=False)

    return df_samples_mags  

def collect_mags_prodigal(s, pckg, working_dir, output_dir):
    """
    Collect Prodigal ORF prediction data for a single sample.
    Copies and renames orf_predicted.faa for each MAG, parses protein stats,
    and joins with MAGs.tsv (checkm) data.
    Returns a DataFrame or None if the sample should be skipped.
    """

    pckg_dir = os.path.join(working_dir, s, pckg)
    file_out = os.path.join(output_dir, f'{s}_MAGs_prodigal.tsv')

    if not os.path.exists(pckg_dir):
        tqdm.write(f'No {s} found at {pckg_dir}. Skipping.')
        return None

    file = os.path.join(pckg_dir, 'MAGs.tsv')
    try:
        df_checkm = pd.read_csv(file, sep='\t')
    except FileNotFoundError:
        tqdm.write(f"Error opening MAGs.tsv: file not found for sample {s}. Skipping.")
        return None

    mags = df_checkm['MAGs'].to_list()
    all_mags_prodigal = []

    for mag in mags:

        n_proteins = 0
        total_protein_length = 0

        faa_file = os.path.join(pckg_dir, mag, 'orf_predicted.faa')
        if not os.path.exists(faa_file):
            tqdm.write(f"orf_predicted.faa not found for MAG {mag} in sample {s}. Skipping.")
            continue

        for record in SeqIO.parse(faa_file, 'fasta'):
            n_proteins += 1
            total_protein_length += len(record.seq)
        mag_file = os.path.dirname(faa_file)

        glab_mag_faa = copy_rename_mag(mag_file, output_dir=output_dir, sample_id=s, glab_signature='GlabCoEvo')

        mag_row = df_checkm[df_checkm['MAGs'] == mag]
        complet = mag_row['Completeness'].values[0]
        contam = mag_row['Contamination'].values[0]

        row = {
            'GLab_faa': glab_mag_faa,
            'MAGs': mag,
            'binID': mag_row['binID'].values[0],
            'sample': s,

            'Marker lineage': mag_row['Marker lineage'].values[0],
            'n_genomes': mag_row['# genomes'].values[0],
            'n_markers': mag_row['# markers'].values[0],
            'n_marker_sets': mag_row['# marker sets'].values[0],
            'Completeness': complet,
            'Contamination': contam,
            'Strain heterogeneity': mag_row['Strain heterogeneity'].values[0],

            'n_proteins': n_proteins,
            'total_protein_length_aa': total_protein_length,
            'mean_protein_length_aa': round(total_protein_length / n_proteins, 2) if n_proteins > 0 else 0,
        }

        all_mags_prodigal.append(row)

    df_samples_mags = pd.DataFrame(all_mags_prodigal)
    #df_samples_mags.to_csv(file_out, sep='\t', index=False)

    return df_samples_mags


def main():

    args = parse_args()

    samples = get_samples(args.samples)
    
    os.makedirs(args.output_dir, exist_ok=True)

    main_file_out_mags = os.path.join(args.output_dir, f'{args.campaigns}_MAGs.tsv')
    main_file_out_mags_prodigal = os.path.join(args.output_dir, f'{args.campaigns}_MAGs_prodigal.tsv')

    all_dataframes_mags = []
    all_dataframes_mags_prodigal = []

    # gathering checkm data for the MAGs
    for s in tqdm(samples, desc="Processing MAGs"):
        df_mag = collect_mags(s, 'mags', args.working_dir, args.output_dir)
        if df_mag is not None:
            all_dataframes_mags.append(df_mag)

    if all_dataframes_mags:
        combined_df_mags = pd.concat(all_dataframes_mags, ignore_index=True)
        combined_df_mags.to_csv(main_file_out_mags, sep='\t', index=False)
        
    # gathering prodigal data for the MAGs
    for s in tqdm(samples, desc="Processing MAGs prodigal"):
        df = collect_mags_prodigal(s, 'mags_prodigal', args.working_dir, args.output_dir)
        if df is not None:
            all_dataframes_mags_prodigal.append(df)

    if all_dataframes_mags_prodigal:
        combined_df_mags_prodigal = pd.concat(all_dataframes_mags_prodigal, ignore_index=True)
        combined_df_mags_prodigal.to_csv(main_file_out_mags_prodigal, sep='\t', index=False)



def parse_args():

    parser = argparse.ArgumentParser("Collect statistics and simple average")
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
    parser.add_argument("-c", "--campaigns",
        help="Name of the campaigns to be used in the output file names.",
        type=str,
        required=True
    )
    parser.add_argument("-o", "--output_dir",
        help="Where the files are wirtten.",
        type=str,
        required=True
    )
    return parser.parse_args()

if __name__ == '__main__':
    main()

# release/
# ├── MAGs/
# │   ├── sample1_mag_1.fa
# │   ├── sample1_mag_2.fa
# │   └── ...
# ├── prodigal/
# │   ├── sample1_mag_1.faa
# │   ├── sample1_mag_2.faa
# │   └── ...
# ├── all_MAGs.tsv          ← combined checkm stats
# └── all_MAGs_prodigal.tsv ← combined prodigal stats

# tar -czf mags_release.tar.gz release/