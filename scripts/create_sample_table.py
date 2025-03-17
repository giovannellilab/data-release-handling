
import os

import pandas as pd

import glob

import argparse

import subprocess

import tempfile


# #  Merge all raw Metagnoems files to a unique folder 
# # Geomosaic requires raw files to be ALL in the same folder
# # Takes files  from /media/edotacca/Loki/sequencing_data/
# # Note: Given the file size, should make a copy to /media/edotacca/Thor/geomosaic_upload
def merge_raw_files(
    dir_samples:str,
    final_dir:str,
    sample_pattern:str
    )->str: 

    pattern_to_sample = os.path.join(dir_samples,sample_pattern)
    all_samples = glob.glob(pattern_to_sample)


    project_name = dir_samples.split(os.path.sep)[-3]
    final_dir = os.path.join(final_dir,project_name)

    FILES = []
    all_rows = []
    pattern_files = '*fq.gz'

    for sample in all_samples:

        files_path = os.path.join(sample,pattern_files)
        all_files = glob.glob(files_path)
        print(all_files)

        relative_files_path = [ "/".join(file.split('/')[-2:]) for file in all_files ]

        FILES.extend(relative_files_path)

        files_names = [ file_p.split('/')[-1] for file_p in all_files ]
        for file in files_names:
            suffix = file.split('_')[-1][0]
            sample_name = file.split('_')[0]

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

    print(FILES)

    frame = pd.concat(all_rows)

    print(frame)

    file_name = f'sample_table_{project_name}.tsv'
    output_file = os.path.join(final_dir,file_name)

    frame.to_csv(output_file, sep = '\t', encoding='utf-8', index=False)
    
    print(f'File {file_name} saved to: {final_dir}')


    # with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmpfile:
    #     # Write absolute paths to the temporary file
    #     tmpfile.write("\n".join(FILES) + "\n")
    #     temp_file_path = tmpfile.name  # Store the temp file path

    # if os.path.exists(final_dir):
    #     print(f'Directory at {final_dir} Found')
    # else:
    #     print(f'Directory not found, creating a new: {final_dir}')
    #     os.makedirs(final_dir)
    
    # # copying file-name list form location to final dir
    # command = ['rsync', '-av', '--no-relative','--files-from='+ temp_file_path, dir_samples, final_dir]
    
    # try:
    #     subprocess.run(command, check=True)
    #     print(f"Files moved to {final_dir}")
    # except subprocess.CalledProcessError as e:
    #     print(f"Error downloading spreadsheet: {e}")

    # os.remove(temp_file_path)


    return output_file



if __name__ == '__main__':


    parser = argparse.ArgumentParser("preprocess_sequences")
    parser.add_argument(
        "-d", "--data_dir",
        help="directory storing sample sequences",
        default=str
    )
    parser.add_argument(
        "-f", "--final_folder",
        help="directory where you wish to distirbute the sequences",
        type=str
    )
    parser.add_argument(
        "-p", "--pattern_sample",
        help="Pattern or prefix for samples",
        type=str
    )

    args = parser.parse_args()

   
    final_folder = merge_raw_files(
        dir_samples = args.data_dir,
        final_dir = args.final_folder,
        sample_pattern = args.pattern_sample
    )

   

