
import os
import pandas as pd
import glob
import argparse
import subprocess
import sys

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

        sample_table = args.sample_table
        if os.path.exists(sample_table):
            print(f'Provided sample table: {os.path.basename(sample_table)} exists')
            df = pd.read_csv(sample_table, sep='	')
            print(df)
        else:
            print(f'Provided sample table or path: {sample_table} does NOT exists')

    else:
        
        print(f'Creating sample table from dir: {args.data_dir}')

        # sample_table = build_sample_table(
        #     dir_samples=args.data_dir,
        #     sample_file=args.subset_samples,
        #     project_name=args.campaign_name,
        #     pattern=args.pattern

        # )
        sample_table = build_sample_table_2(
            dir_samples=args.data_dir,
            sample_file=args.subset_samples,
            project_name=args.campaign_name,
            pattern=args.pattern

        )
    ibisco_uploader(
        sample_table=sample_table,
        dir_samples=args.data_dir,
        project_name=args.campaign_name,
        user_name=args.ibisco_user,
        dry_run=args.dry_run
    )

def parse_args():

    parser = argparse.ArgumentParser("Upload raw seqs AND sample_table.tsv to IBISCO server")
    parser.add_argument(
        "-d", "--data_dir",
        help="directory storing sample sequences",
        default=str
    )
    parser.add_argument(
        "-t", "--sample_table",
        help="provides geomosaic sample table",
        type=str,
        default=None
    )
    parser.add_argument(
        "-c", "--campaign_name",
        help="Provide campaign name: Ex ARG23 or TEST",
        type=str
    )
    parser.add_argument(
        "-p", "--pattern",
        help="Sequence pattern to search, default: *fq.gz",
        type=str,
        default='*fq.gz'
    )
    parser.add_argument(
        "-s", "--subset_samples",
        help="provide a text file with each line a sample",
        type=str,
        default=None
    )
    parser.add_argument(
        "-u", "--ibisco_user",
        help="Provide user account",
        type=str,
        default='dgiovannelli'
    )
    parser.add_argument(
        "-z","--dry_run",
        help="Will attempt a dry run without uploading files",
        action="store_true"
    )

    return parser.parse_args()


def subset(
        sample_file:str,
)-> list:
    """
    Select samples specified in a text file
    Args:
        textfile (str): sample x row.
    Returns:
        list: samples
    """
    samples = []

    with open(sample_file,'r')as reader:
        for line in reader:
            sample = line.strip()
            samples.append(sample)
    
    return samples


def build_sample_table_2(dir_samples: str, sample_file: str, project_name: str, pattern: str) -> str: 

    if not pattern.startswith('*'):
        search_pattern = f"**/*{pattern}"
    else:
        search_pattern = f"**/{pattern}"

    print(f"Searching for: {search_pattern} in {dir_samples}")

    all_files_path = glob.glob(os.path.join(dir_samples, search_pattern), recursive=True)
    all_files_path = [f for f in all_files_path if os.path.isfile(f)]

    if not all_files_path:
        print(f"ERROR: No files found matching {pattern} in {dir_samples}")
        sys.exit(1)

    # 3. Group files by Sample Name
    samples_dict = {}

    for f_path in sorted(all_files_path):

        fname = os.path.basename(f_path)
                parts = fname.split('_')
        if len(parts) < 2:
            continue
            
        # Example: HTB_S21 (from HTB_S21_R1)
        sample_id = "_".join(parts[:-1])
        # Example: R1.fastq.gz
        suffix = parts[-1]

        if sample_id not in samples_dict:
            samples_dict[sample_id] = {'r1': None, 'r2': None}

        if 'R1' in suffix:
            samples_dict[sample_id]['r1'] = fname
        elif 'R2' in suffix:
            samples_dict[sample_id]['r2'] = fname
        else:
            print(f"WARNING: File {fname} does not match expected R1/R2 pattern at the end, skipping.")
            continue

    all_rows = []
    for s_id, reads in samples_dict.items():
        row = pd.DataFrame([{
            'r1': reads['r1'],
            'r2': reads['r2'],
            'sample': s_id
        }])
        all_rows.append(row)

    frame = pd.concat(all_rows, ignore_index=True)
    print("\nGenerated Sample Table:")
    print(frame)

    file_name_table = f'sample_table_{project_name}.tsv'
    sample_table_file = os.path.join(dir_samples, file_name_table)
    frame.to_csv(sample_table_file, sep='\t', index=False)
    print(f'STEP[2] SAVING File {file_name_table} to: {dir_samples}')

    
    return sample_table_file



def ibisco_uploader(
    sample_table: str,
    dir_samples: str,
    project_name: str,
    user_name: str,
    dry_run: bool = False  # default, used if not overridden
    ):
    # # this function upload files tored in the sample_table to our SERVER # #
    # # the files are copied in the assgined 'final_dir' # #
    sample_df = pd.read_csv(sample_table, sep = '\t')
    if not os.path.isdir(dir_samples):
        print(f'Path not valid: {dir_samples}')
        sys.exi(1)

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
    
    HOST_name = 'ibiscohpc-ui.scope.unina.it'
    user_server = '@'.join([user_name,HOST_name])

    if user_name == 'dgiovannelli':
        remote_path = f'/ibiscostorage/GiovannelliLab/raw/{project_name}'
    else:
        remote_path = f'/ibiscostorage/{user_name}/raw/{project_name}'

    remote_host = f'{user_server}:{remote_path}'
    
    pre_command = ['ssh', user_server, 'mkdir', '-p', remote_path]

    if dry_run:
        print('NOT uploadnig -> dry-run activated \n Samples table built')
        command = [
                    'rsync',
                    '--dry-run',
                    '-av', 
                    '--progress',
                    '--no-relative',
                    '--files-from=-',
                    '/',
                    remote_host
                    ]
    else:
        command = [
                    'rsync', 
                    '-av', 
                    '--progress',
                    '--no-relative',
                    '--files-from=-',
                    '/', 
                    remote_host
                    ]
    try:
        if dry_run:
            subprocess.run(pre_command, check=True)
            subprocess.run(command, input=files_payload, text=True, check=True)
        else:
            print(f'Data will be uploaded to: {remote_host}')
            subprocess.run(pre_command, check=True)
            subprocess.run(command, input=files_payload, text=True, check=True)
            print(f'\n STEP[3] Files uploaded to: {remote_host}')
    
    except subprocess.CalledProcessError as e:
        print(f"\n Error while copying file: {e}")



if __name__ == '__main__':

    main()