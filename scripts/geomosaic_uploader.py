
import os

import pandas as pd

import glob

import argparse

import subprocess

import tempfile

import sys


# NOTE:
# to create the text file with all dirs names, i will recommend to place yoursefl in the parent dir
# storing all the sample's sub dirs and run the follwing command:
# find . -mindepth 1 -maxdepth 1 -type d -exec basename {} \; > dirs.txt


# #  Merge all raw Metagnoems files to a unique folder 
# # Geomosaic requires raw files to be ALL in the same folder
# # Takes files  from /media/edotacca/Loki/sequencing_data/



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



def build_sample_table(
    dir_samples:str,
    sample_file:str
    )->str: 

    #pattern_to_sample = os.path.join(dir_samples,sample_pattern)
    pattern_files = '*fq.gz'

    if sample_file:

        samples = subset(
            sample_file=sample_file
            )
        all_samples = [glob.glob(os.path.join(dir_samples,file)) for file in samples]

    else:
        all_dir_samples = os.path.join(f'{dir_samples}/**/')
        all_samples = [dire for dire in glob.glob(all_dir_samples)]
        print(all_samples)

    project_name = dir_samples.split(os.path.sep)[-2]
    FILES = []
    all_rows = []
    ## WGs file extension

    for sample in all_samples:
    
        
        sample = sample[0]
        print('sample:',sample)
        #selecting samples's files
        files_path = os.path.join(sample,pattern_files)
        all_files = glob.glob(files_path)

        #extracting sample/file path
        relative_files_path = [ "/".join(file.split('/')[-2:]) for file in all_files ]
        FILES.extend(relative_files_path)

        # extracting file names
        files_names = [ file_p.split('/')[-1] for file_p in all_files ]
        print(files_names)

        for file in files_names:
            #suffix could be changed 
            suffix = file.split('_')[-1][0]
            # if sample_name is cannonical G100; jsut splice with [0]
            sample_name = "_".join(file.split('_')[:3])
            #sample_name = file.split('_')[0]

            if suffix == '1':
                file_forward = file
            if suffix == '2':
                file_reverse = file

        row = pd.Series({
                'r1' : file_forward ,
                'r2' : file_reverse ,
                'sample' : sample_name
        }).to_frame().T
        all_rows.append(row)

    frame = pd.concat(all_rows)
    print(frame)

    # # saving sample table file

    file_name_table = f'sample_table_{project_name}.tsv'
    sample_table_file = os.path.join(dir_samples,file_name_table)
    path_to_sample_table = sample_table_file.split('/')[-1]

    frame.to_csv(sample_table_file, sep = '\t', encoding='utf-8', index=False)
    print(f'STEP[2] SAVING File {file_name_table} to: {dir_samples}')

    # # saving sample file paths

    file_name_paths = f'file_paths_{project_name}.txt'
    file_paths = os.path.join(dir_samples,file_name_paths)
    # wriring to a temporary file for rsync upload
    with open(file_paths,'w') as writer:
        writer.write('\n'.join(f'/{path}' for path in FILES) + '\n')
        writer.write(f'{path_to_sample_table}')

    
    return file_name_table, file_paths


def ibisco_uploader(
    sample_table: str,
    file_paths: str,
    dir_samples: str,
    user_name: str,
    dry_run: bool = False  # default, used if not overridden
    ):
    # # this function upload files tored in the sample_table to our cluster # #
    # # the files are copied in the assgined 'final_dir' # #


    # # reading table
    sample_table_name = os.path.basename(sample_table)
    project_name = os.path.dirname(dir_samples).split('/')[-1]
    print('project_name',project_name)
    
    ## replace the user anme with your own account
    HOST_name = 'ibiscohpc-ui.scope.unina.it'
    user_server = '@'.join([user_name,HOST_name])

    if user_name == 'dgiovannelli':
        remote_path = f'GiovannelliLab/raw/{project_name}'
    else:
        remote_path = f'{user_name}/raw/{project_name}'

    final_dir = f'{user_server}:/ibiscostorage/{remote_path}'
    
    print(f'Data will be uploaded to: {final_dir}')
   
    if dry_run:

        command = [
                    'rsync',
                    '--dry-run',
                    '-av', 
                    '--progress',
                    '--no-relative',
                    '--files-from='+ file_paths, 
                    dir_samples, 
                    final_dir
                    ]

    else:
        command = [
                    'rsync', 
                    '-av', 
                    '--progress',
                    '--no-relative',
                    '--files-from='+ file_paths, 
                    dir_samples, 
                    final_dir
                    ]

    try:
        subprocess.run(command, check=True)
        print(f'\n STEP[3] Files uploaded to: {final_dir}')
        os.remove(file_paths)
    
    except subprocess.CalledProcessError as e:
        print(f"\n Error while copying file: {e}")




if __name__ == '__main__':


    parser = argparse.ArgumentParser("preprocess_sequences")
    parser.add_argument(
        "-d", "--data_dir",
        help="directory storing sample sequences",
        default=str
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
        "--dry-run",
        help="Will attempt a dry run without uploading files",
        action="store_true"
    )

    args = parser.parse_args()


    sample_table, paths_file = build_sample_table(
        dir_samples=args.data_dir,
        sample_file=args.subset_samples,
    )
    ibisco_uploader(
        sample_table=sample_table,
        file_paths=paths_file,
        dir_samples=args.data_dir,
        user_name=args.ibisco_user     
    )
   

