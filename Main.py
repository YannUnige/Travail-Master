import Profiles as Pr
import Plotting as Pl
import os
import Matrix as Mat
import sys
import Hits_List as HL
import shutil
import pandas as pd
import glob
import argparse
import time
# Sorts the input group by evalue in ascending order and drops duplicates based on 'id', keeping the first occurrence.
def process_group(group):
        return group.sort_values(by='evalue', ascending=True).drop_duplicates(subset='id', keep='first')

#  Remove duplicated protein data from multiple CSV files and save the processed data back to each file.
def Remove_duplicated_protein(Results_path):
    pattern = f"{Results_path}/**/*_Filtered_Sequence_*.csv"
    files = glob.glob(pattern, recursive=True)

    # Prepare a list to hold dataframes
    dfs = []
    for file in files:
        df = pd.read_csv(file)
        # Extract model_name 
        model_name = file.split('/')[-2] 
        df['model'] = model_name
        dfs.append(df)
    
    # Combine all DataFrames 
    combined_df = pd.concat(dfs)
    #print("combined_df\n",combined_df.to_string())
    mean_evalues = combined_df.groupby(['id', 'model'])['evalue'].mean().reset_index()
    #print("mean_evalues\n", mean_evalues.to_string())

    best_model_per_id = mean_evalues.loc[mean_evalues.groupby('id')['evalue'].idxmin()]
    #print("-------best_model_per_id\n",best_model_per_id.to_string())
    best_model_per_id.rename(columns={'evalue': 'mean_evalue'}, inplace=True)

    processed_df = pd.merge(combined_df, best_model_per_id, on=['id', 'model'], how='inner')

    #print("processed_df\n",processed_df.to_string())


    # Separate and save the processed data
    for file in files:
        model_name = file.split('/')[-2]  
        model_df = processed_df[processed_df['model'] == model_name].copy()
        model_df.drop(columns=['model'], inplace=True)  # Clean up by removing the model column
        model_df.drop(columns=['mean_evalue'], inplace=True)
        model_df.to_csv(file, index=False)


#Remove model with non significant hits from the specified parent directory by deleting subdirectories that do not contain the target folder.
def Remove_non_hits(parent_directory, target_folder):
    # Iterate over each item in the parent directory
    for item in os.listdir(parent_directory):
        item_path = os.path.join(parent_directory, item)
        if os.path.isdir(item_path):
            target_folder_path = os.path.join(item_path, target_folder)
            # Check if the "Hits_Profiles_Filtered" folder exists within the subdirectory
            if not os.path.exists(target_folder_path):
                # If not found, delete the subdirectory
                shutil.rmtree(item_path)
                
# Parse command-line arguments         
def parse_arguments():
    parser = argparse.ArgumentParser(description="Process proteome data.")
    parser.add_argument("Proteome_path", type=str, help="Path to the proteome directory.")
    parser.add_argument("Model_path", type=str, help="Path to the directory containing the model.")
    parser.add_argument("Bank", type=str, help="Bank to use for analysis.")
    parser.add_argument("--evalue", type=float, default=0.012, help="E-value threshold for hits (default: 0.012).")    
    
    # Check the number of arguments provided
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    
    return parser.parse_args()
# Extracts and returns the organism name from the proteome path
def get_organism_name(path):
    path_segments = path.split('/')
    return path_segments[2] if len(path_segments) >= 3 else "Unknown"


if __name__ == "__main__":
    start = time.time()
    args = parse_arguments()

    Proteome_path = args.Proteome_path
    Model_path = args.Model_path
    Bank = args.Bank
    Evalue = args.evalue

    organism_name = get_organism_name(args.Proteome_path)

    print("Proteome Path:", args.Proteome_path)
    print("Directory Model:", args.Model_path)
    print("Bank:", args.Bank)
    print("Organism Name:", organism_name)
 

    Result_Path = f"Results_{organism_name}_{Bank}_{Evalue}"
    # Checking data folders
    if not os.path.exists(Proteome_path):
        print("Proteome path not found")
        exit()
    if not os.path.exists(Model_path):
        print("Model path not found")
        exit()
    
    # Checking results folder
    if os.path.exists(Result_Path):
        print("Results folder already exists, removing it ? (y/n)")
        answer = input()
        if answer.lower() == 'y':
            shutil.rmtree(Result_Path)
        else:
            exit()
    
    # Iterate over Model to apply hmmer search
    for model_filename in os.listdir(Model_path):
        
        Model = os.path.join(Model_path, model_filename)
        
        print("Start of the search for",model_filename)    
        Model_name = Model.split("/")[-1].replace(".hmm","")

        Results_path_Model = f"{Result_Path}/{Model_name}"
        print("Hmmer Search")
        Filtered_sequences = Pr.profiles(Proteome_path,Model,Bank,Evalue,Results_path_Model,Model_name)
    
    Remove_duplicated_protein(Result_Path)

    # Plotting for each model
    for model_filename in os.listdir(Model_path):
        
        Model = os.path.join(Model_path, model_filename)

        Model_name = Model.split("/")[-1].replace(".hmm","")
        print(Model_name)

        Results_path_Model = f"{Result_Path}/{Model_name}"
        Filtered_sequences = f"{Results_path_Model}/{Model_name}_Filtered_Sequence_{Bank}.csv"
        print("HeatMap generation")
        if not os.path.exists(Filtered_sequences):
            print("Filtered sequences not found")
            exit()
        Pl.Plot(Results_path_Model,Filtered_sequences,Model_name,Bank)
        print("End of the search for",model_filename)
        print("\n") 
   
 
    HL.List_Hits_txt(Result_Path,organism_name)
    Mat.Matrix_Hits(Proteome_path,Model_path,Result_Path,Bank,organism_name)
    Remove_non_hits(Result_Path,"Hits_Profiles_Filtered")
    end = time.time()
   
        
    with open("Time.txt", "a") as f:
        f.write(f"Total time taken for {organism_name} {Bank}: {end - start:.2f} seconds \n")


    
    