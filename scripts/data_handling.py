import subprocess

import pandas as pd

import argparse

import os

import glob

import pickle

from datetime import datetime

import pprint


def check_sanity()-> None:

    # # check for the presence of spreadsheet

    # # check for the presence of final directories

    return None


def download_spreadsheet(
        spreadsheet_id:str, 
        sheet_id:str
        )-> str:
    
    output_file = "sample_database_biosample_id.csv"
    template_dir = 'template'
    output_wget= os.path.join(template_dir,output_file)

    url = f"https://docs.google.com/spreadsheets/d/{spreadsheet_id}/gviz/tq?tqx=out:csv&gid={sheet_id}"
    command = ["wget","-N" "-O", output_wget, url]
    
    # try:
    #     subprocess.run(command, check=True)
    #     print(f"Spreadsheet downloaded successfully (if updated) as {output_file}")
    # except subprocess.CalledProcessError as e:
    #     print(f"Error downloading spreadsheet: {e}")

    return output_wget



def merger_folder(
        data_release: str
) -> str:

    # # grab data release dir
    release = os.path.basename(data_release)
    source_directory = os.path.join(data_release,f'result_{release}')
    print(source_directory)
    clean_data = os.path.join(source_directory,'00.CleanData')
    raw_data = os.path.join(source_directory,'01.RawData')

    # create a combined directry for clean and raw data
    combined = os.path.join(source_directory,'02.Combined')

    # Ensure the combined directory exists
    os.makedirs(combined, exist_ok=True)

    #Use rsync to merge directories (preserving structure & skipping duplicates)
    try:
        subprocess.run(['rsync', '-a', '--ignore-existing', raw_data + '/', combined], check=True)
        subprocess.run(['rsync', '-a', '--ignore-existing', clean_data + '/', combined], check=True)
        print(f"Successfully merged folders into {combined}")
    except subprocess.CalledProcessError as e:
        print(f"Error during rsync: {e}")


    return combined

def map_samples(
        data_release:str,
        csv_file:str

):  
    # # grab data release dir

    # release = os.path.basename(data_release)
    # source_directory = os.path.join(data_release,f'result_{release}')
    # print(source_directory)
    samples_dir = os.path.join(data_release,
                               '01.RawData'
                               )

    dataframe = pd.read_csv(csv_file, 
                            index_col=False)
    #print(dataframe)

    fields = ['ExpID','Sample_name',
              'amplicon_Univ V45 (U)']

    dataframe = dataframe\
        .dropna(subset=[fields[2]])

    data_df = dataframe[fields]
    associations = {}

    for i,row in data_df.iterrows():
        if row.iloc[2] not in associations.keys():
            associations[row.iloc[2]] = row.iloc[0]

    if not os.path.exists(samples_dir):
        print(f'Error {samples_dir} not correctly inputed')
        return 

    sample_map = {}
    pattern = samples_dir + '/G*'

    for sample_dir in glob.glob(pattern):
        
        sample = os.path.basename(sample_dir)

        if sample in associations.keys():
            key = associations[sample]

            if key not in sample_map:
                sample_map[key] = [sample]
            elif key in sample_map:
                sample_map[key].append(sample)

    # pd.Series({
    #         'Batch_ID':,
    #         'download_date':,
    # }).frame().T

    print('Campaigns:',sample_map.keys())
    print('Samples:',len(samples_dir))
    #print(sample_map.items())
    for key,value in sample_map.items():
        print(f'> {key} --> {value}')

    # OPen the csv file as a dataframe
    # check for the rpesence of same Batch_ID,
    # if TRue exit
    # else -> append new row 
    ### return a csv file with info:
    # - batch_ID
    # - date
    # - all samples

    return sample_map



def save_release_info(
        map_sample: dict, 
        experiment_type: str,
        data_release: str
        ):
    now = datetime.now()
    #formatted_date = now.strftime("%d-%m-%Y")
    #formatted_date = '12-02-2025'
    formatted_date = '21-02-2025'
    # info = { formatted_date: map_sample.values()}
    info = {}
    info[formatted_date] = map_sample.values

    template_dir = "template"
    pickle_file = "data_release_2025.pickle"
    location_file = os.path.join(template_dir,pickle_file)

    print(info)

    # Check if pickle file exists
    if os.path.exists(location_file):
        try:
            with open(location_file, "rb") as handle:
                existing_data = pickle.load(handle)  # Load existing data

        except (EOFError, pickle.UnpicklingError):
            existing_data = {}  # If the file is corrupted or empty, initialize with empty dict
    else:
        existing_data = {}

    # Merge new data with existing data
    # existing_data.update(info)

    # # Save updated dictionary back to pickle file
    # with open(pickle_file, "wb") as handle:
    #     pickle.dump(existing_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Data saved successfully to {location_file}")



def distribute_samples(
        final_folder:str,
        data_release:str,
        map_samples:dict,
        experiment_type:str
):
    final_fo = final_folder + '/*'
    campaigns = glob.glob(final_fo)
    print(campaigns)
    release = 'X204SC25036218-Z01-F001'

    # # source dir
    data_release_dir = os.path.join(
        data_release,
        '02.Combined'
        )
    report_path = os.path.join(data_release,f'report_{release}')

    for campaign,samples in map_samples.items():
        #check if campaign exists in final destination:
        current_dir = os.path.join(final_folder,campaign)
        print(f'Current directory: {current_dir}')
        #subprocess.call('ls', shell=True, cwd=current_dir)

        if not os.path.exists(current_dir):
            print(f'Error {current_dir} not correctly inputed')
            return 
        
        else:
            subdir = os.path.join(current_dir,experiment_type)
            # cehck if dir already exists
            if not os.path.exists(subdir):
                os.makedirs(subdir)
                print(f'{subdir} created')
            else:
                print(f'{subdir} already existed')
        #retrieve all samples associate to the current campaigns
        new_samples = [os.path.join(data_release_dir,sample) for sample in samples ]

        for s in new_samples:
            print(s)

        try:
            subprocess.run(['rsync',
                            '-av', 
                            '--ignore-existing'] 
                            + new_samples 
                            + [subdir], 
                            check=True)

            subprocess.run(['rsync',
                            '-av', 
                            '--ignore-existing']
                            + [report_path]
                            + [subdir],
                            check=True)
        
        except subprocess.CalledProcessError as e:
            print(f"Error during rsync: {e}")            

        print(f'{os.path.basename(current_dir)} --> copied --------------------')

if __name__ == '__main__':

    parser = argparse.ArgumentParser("preprocess_sequences")
    parser.add_argument(
        "-e", "--experiment_type",
        help="type of files: 16_S or 18_S or ITS",
        type=str
    )
    parser.add_argument(
        "-d", "--data_release",
        help="directory storing new sample sequences",
        default=str
    )
    parser.add_argument(
        "-f", "--final_folder",
        help="directory where you wish to distirbute the new sequences",
        type=str
    )
    args = parser.parse_args()

    
    # Replace with actual values
    SPREADSHEET_ID = "1ecIerj5tziIqG374FzAQyCjSdhwaxxJ2Q2H4ZdJuLjc"
    SHEET_ID = "1748512334"
    OUTPUT_FILE = "sample_database_biosample_id.csv"


    output_file = download_spreadsheet(
            spreadsheet_id =SPREADSHEET_ID,
            sheet_id=SHEET_ID,
                        )
    
    mapped = map_samples(
        data_release=args.data_release,
        csv_file=output_file
    )

    # combined_dir = merger_folder(
    #     data_release=args.data_release
    # )

    # save_release_info(
    #     map_sample = mapped,
    #     experiment_type=args.experiment_type,
    #     data_release=args.data_release
    #     )
    
    distribute_samples(
        final_folder =args.final_folder,
        data_release=args.data_release,
        map_samples=mapped,
        experiment_type=args.experiment_type
        )
    