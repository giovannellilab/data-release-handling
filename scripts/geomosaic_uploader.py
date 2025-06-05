
import os

import pandas as pd

import glob

import argparse

import subprocess

import tempfile


# python3 create_sample_table.py  
# -d /media/edotacca/Loki/sequencing_data/ISL22/Metagenomes/ 
# -f /dirs.txt


# NOTE:
# to create the text file with all dirs names, i will recommend to place yoursefl in the parent dir
# storing all the sample's sub dirs and run the follwing command:
# find . -mindepth 1 -maxdepth 1 -type d -exec basename {} \; > dirs.txt


# #  Merge all raw Metagnoems files to a unique folder 
# # Geomosaic requires raw files to be ALL in the same folder
# # Takes files  from /media/edotacca/Loki/sequencing_data/


def subset(
        sample_file:str,
):
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


    if samples:
        all_samples = [glob.glob(os.path.join(dir_samples,file)) for file in samples]
    else:
        all_samples = glob.glob(dir_samples)


    project_name = dir_samples.split(os.path.sep)[-3]
    #final_dir = os.path.join(final_dir,project_name)

    FILES = []
    all_rows = []
    ## WGs file extension
    pattern_files = '*fq.gz'

    for sample in all_samples:
        
        sample = sample[0]
        #selecting samples's files
        files_path = os.path.join(sample,pattern_files)
        all_files = glob.glob(files_path)

        #extracting sample/file path
        relative_files_path = [ "/".join(file.split('/')[-2:]) for file in all_files ]
        FILES.extend(relative_files_path)

        # extracting file names
        files_names = [ file_p.split('/')[-1] for file_p in all_files ]

        for file in files_names:
            suffix = file.split('_')[-1][0]
            # if sample_name is cannonical G100; jsut splice with [0]
            sample_name = "_".join(file.split('_')[:3])

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

    # # saving file

    file_name_table = f'sample_table_{project_name}.tsv'
    sample_table_file = os.path.join(dir_samples,file_name_table)

    frame.to_csv(file_name_table, sep = '\t', encoding='utf-8', index=False)
    print(f'STEP[2] SAVING File {file_name_table} to: {dir_samples}')


    file_name_paths = f'file_paths_{project_name}.txt'
    file_paths = os.path.join(dir_samples,file_name_paths)
    # wriring to a temporary file for rsync upload
    with open(file_paths,'w') as writer:
        writer.write('\n'.join(f'/{path}' for path in FILES) + '\n')
        writer.write('\n'.join(f'{sample_table_file}'))
    
    return file_name_table, file_paths


def ibisco_uploader(
    sample_table : str,
    file_paths : str,
    dir_samples: str,
    user_name : str,
    ask_upload: bool = False,
    ):
    # # this function upload files tored in the sample_table to our cluster # #
    # # the files are copied in the assgined 'final_dir' # #


    # # reading table
    sample_table_name = os.path.basename(sample_table)    
    project_name = os.path.dirname(dir_samples).split('/')[-2]

    
    ## replace the user anme with your own account
    HOST_name = 'ibiscohpc-ui.scope.unina.it'
    user_server = '@'.join([user_name,HOST_name])

    final_dir = f'{user_server}:/ibiscostorage/GiovannelliLab/raw/{project_name}'

    # # # ask
    # if ask_upload:
    #     reply = input("Proceed with uploading files? Type 'yes' or 'no': ").strip().lower()
    #     if reply != 'yes':
    #         print("Upload aborted by user.")
    #         sys.exit(1)  # Exits with error code 1

    # copying file-name list form location to final dir
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
    except subprocess.CalledProcessError as e:
        print(f"\n Error while copying file: {e}")
    os.remove(file_paths)




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
        "--confirm-upload",
        help="Ask user whether to proceed with file upload (yes/no)",
        action="store_true"
    )
    args = parser.parse_args()

    samples = subset(
        sample_file=args.subset_samples
        )
 # parser.add_argument(
    #     "-p", "--pattern_sample",
    #     help="Pattern or prefix for samples",
    #     type=str,
    #     default='G*'
    # )
    sample_table, paths_file = build_sample_table(
        dir_samples = args.data_dir,
        sample_file=samples,
    )
    ibisco_uploader(
        sample_table=sample_table,
        file_paths=paths_file,
        dir_samples=args.data_dir,
        ask_upload=args.confirm_upload,
        user_name=args.ibisco_user     
    )
   

