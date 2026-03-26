
import os
import pandas as pd
import glob
import argparse
import subprocess
import sys
import re
from pathlib import Path

### AUTHOR: EDOARDO TACCALITI ###

# This script uses a sample table (the one used by Geomosaic setup ) to upload WGS seqeunces
# in the raw/ partition for initial setup.
# In the case NO smape table is passed, the script tries to generate one from the path passed via
# command line parameters

# ASK for HELP:
# conda activate ena
# cd data-release-handling
# python scripts/geomosaic_uploader.py -h
 
def main():

    args = parse_args()

    if args.sample_table:

        sample_table_file = args.sample_table

        if os.path.exists(sample_table_file):
            print(f'Provided sample table: {os.path.basename(sample_table_file)} exists')
            sample_table = table_checks(sample_table_file)
        else:
            print(f'Provided sample table or path: {sample_table} does NOT exists')

    else:
        
        print(f'Creating sample table from dir: {args.data_dir}')
        sample_table = build_sample_table(
            dir_samples=args.data_dir,
            project_name=args.campaign_name

        )
    ibisco_uploader(
        sample_table=sample_table,
        dir_samples=args.data_dir,
        project_name=args.campaign_name,
        remote_host=args.remote_host,
        remote_path=args.remote_path,
        dry_run=args.dry_run
    )

def parse_args():

    parser = argparse.ArgumentParser("Upload raw seqs AND sample_table.tsv to IBISCO server")
    parser.add_argument(
        "-d", "--data_dir",
        help="Provide directory ABSOLUTE path storing sample sequences",
        default=str
    )
    parser.add_argument(
        "-t", "--sample_table",
        help="Provide geomosaic ABSOLUTE path to sample table",
        type=str,
        default=None
    )
    parser.add_argument(
        "-H", "--remote_host",
        help="Provide user & remote server hostname or IP (e.g. dgiovannelli@ibiscohpc-ui.scope.unina.it)",
        type=str,
        default='ibisco',
        required=True
    )
    parser.add_argument(
        "-p", "--remote_path",
        help="Provide absolute remote path on the server (e.g. /ibiscostorage/GiovannelliLab/raw)",
        type=str,
        required=True
    )
    parser.add_argument(
        "-c", "--campaign_name",
        help="Provide output folder name to store files in remote: Ex ARG23 or TEST",
        type=str,
        required=True
    )
    parser.add_argument(
        "-z","--dry_run",
        help="Will attempt a dry run without uploading files",
        action="store_true"
    )

    return parser.parse_args()



def build_sample_table(dir_samples: str, project_name: str) -> str: 

    base_path = Path(dir_samples).resolve()

    if not base_path.is_dir():
        sys.exit(f"ERROR: Directory '{base_path}' does not exist or is not a directory.")

    pattern = r'(\.fastq|\.fq)(\.gz)?$'
    regex_pattern = re.compile(pattern, re.IGNORECASE)

    all_files = [
        f for f in base_path.rglob("*")
        if f.is_file() and regex_pattern.search(f.name)
    ]

    print(f"--- Searching recursively in {base_path} for: {pattern} ---")

    if not all_files:
        sys.exit(f"ERROR: No files found in {base_path} matching {pattern}")

    samples_dict = {}
    
    # Common patterns for R1 and R2
    # This covers _R1, _1, .1, _R1_001, etc.
    r1_regex = re.compile(r'(_R1(_\d+)?|_1)(?=\.fastq|\.fq)(\.fastq|\.fq)(\.gz)?$', re.IGNORECASE)
    r2_regex = re.compile(r'(_R2(_\d+)?|_2)(?=\.fastq|\.fq)(\.fastq|\.fq)(\.gz)?$', re.IGNORECASE)

    for f_path in sorted(all_files):
        fname = f_path.name
        
        if r1_regex.search(fname):
            sample_id = r1_regex.sub('', fname)
            direction = 'r1'
        elif r2_regex.search(fname):
            sample_id = r2_regex.sub('', fname)
            direction = 'r2'
        else:
            # Fallback for unpaired or weirdly named files
            sample_id = fname.split('.')[0]
            direction = 'r1' 

        if sample_id not in samples_dict:
            samples_dict[sample_id] = {'r1': None, 'r2': None}
        
        samples_dict[sample_id][direction] = fname

    rows = []
    for s_id, reads in samples_dict.items():
        rows.append({'r1': reads['r1'], 'r2': reads['r2'], 'sample': s_id})
    
    df = pd.DataFrame(rows).sort_values('sample')
    
    output_path = base_path / f'sample_table_{project_name}.tsv'
    df.to_csv(output_path, sep='\t', index=False)
    
    print("\n--- Generated Sample Table ---\n")
    print(df.to_string(index=False))
    
    return str(output_path)


def table_checks(filename):

    _, extension = os.path.splitext(filename)
    file_format = extension.lower().replace(".", "")

    if file_format == "tsv":
        rawdf = pd.read_csv(filename, sep="\t")
    elif file_format == "csv":
        rawdf = pd.read_csv(filename, sep=",")
    else:
        rawdf = pd.read_excel(filename)
    
    assert "r1" in list(rawdf.columns), f"\n\n ERROR in Column 'r1' not present in the provided table. It should contains three columns, with the following header (all lower-case): r1 r2 sample"
    assert "r2" in list(rawdf.columns), f"\n\n ERROR in Column 'r2' not present in the provided table. It  contains three columns, with the following header (all lower-case): r1 r2 sample"
    assert "sample" in list(rawdf.columns), f"\n\n ERROR in Column 'sample' not present in the provided table. It should contains three columns, with the following header (all lower-case): r1 r2 sample"
    
    print('\n',rawdf.to_string(index=False),'\n')

    return filename


def ibisco_uploader(
    sample_table: str,
    dir_samples: str,
    project_name: str,
    remote_host: str,
    remote_path: str,
    dry_run: bool = False  # default, used if not overridden
    ):
    # # this function upload files tored in the sample_table to our SERVER # #
    # # the files are copied in the assgined 'final_dir' # #
    sample_df = pd.read_csv(sample_table, sep = '\t')
    
    base_path = Path(dir_samples).resolve()

    if not base_path.is_dir():
        sys.exit(f"ERROR: Directory '{base_path}' does not exist or is not a directory.")

    list_all_paths = []
    for row in sample_df.itertuples():
        name_r1 = str(row.r1).strip()
        name_r2 = str(row.r2).strip()

        r1_path = os.path.join(dir_samples, f"**/{name_r1}")
        r2_path = os.path.join(dir_samples, f"**/{name_r2}")
        
        file_1 = glob.glob(r1_path, recursive=True)
        file_2 = glob.glob(r2_path, recursive=True)

        list_all_paths.extend(file_1 + file_2)

    list_all_paths.extend([sample_table])
    files_payload = "\n".join(list_all_paths)
    
    #HOST_name = 'ibiscohpc-ui.scope.unina.it'
    #user_server = '@'.join([user_name,HOST_name])

    full_remote_path = f'{remote_path}/{project_name}'
    full_remote_host = f'{remote_host}:{full_remote_path}'

    pre_command = ['ssh', remote_host, 'mkdir', '-p', full_remote_path]

    if dry_run:
        print('NOT uploadnig -> dry-run activated')
        command = [
                    'rsync',
                    '--dry-run',
                    '-av', 
                    '--progress',
                    '--no-relative',
                    '--files-from=-',
                    '/',
                    full_remote_host
                    ]
    else:
        command = [
                    'rsync', 
                    '-av', 
                    '--progress',
                    '--no-relative',
                    '--files-from=-',
                    '/', 
                    full_remote_host
                    ]
    try:
        if dry_run:
            print(f'Data will be uploaded to: {full_remote_host}')
            subprocess.run(pre_command, check=True)
            subprocess.run(command, input=files_payload, text=True, check=True)
            print(f'\n STEP[3] Attempted Dry-run to: {remote_host}')

        else:
            print(f'Data is being uploaded to: {full_remote_host}')
            subprocess.run(pre_command, check=True)
            subprocess.run(command, input=files_payload, text=True, check=True)
            print(f'\n STEP[3] Files uploaded to: {remote_host}')
    
    except subprocess.CalledProcessError as e:
        print(f"\n Error while copying file: {e}")



if __name__ == '__main__':

    main()