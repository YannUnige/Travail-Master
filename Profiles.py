import os
import subprocess
import shutil
import pandas as pd
from Bio import SearchIO, SeqIO
import time
from multiprocessing import Pool
from multiprocessing import cpu_count


attributs_hits = ["bitscore", "evalue", "id","sequence","assembly","genomic_location","description"] # Attributes extracted

#HMMER search in all proteomes
def Hmmer_Search(model_file, data_genome,Hmmer_output_path,Bank,Model_name):
    if not os.path.exists(Hmmer_output_path):
        os.makedirs(Hmmer_output_path)

    for root, dirs, files in os.walk(data_genome):
        for dir_name in dirs:
            if dir_name.startswith(Bank):
                database_path = os.path.join(root, dir_name, "protein.faa")
                if os.path.exists(database_path):
                    result_file = f"{Hmmer_output_path}/{Model_name}_{dir_name}_Profile.hmmer"
                    with open(result_file, 'w') as output_file:
                        command = ["hmmsearch","--notextw","--domtblout",result_file, model_file, database_path]
                        subprocess.run(command, check=True, stdout=output_file)
    return




#Filter the HMMER hits
def Hmmer_Hits_Filter(Input_folder, Ouput_folder):
    if not os.path.exists(Ouput_folder):
        os.makedirs(Ouput_folder)

    for root, dirs, files in os.walk(Input_folder):
        for file_name in files:
            result_file = os.path.join(root, file_name)
            with open(result_file, 'r') as f:
                lines = f.readlines()
                for line in lines :
                    if "target name" in line:
                        shutil.copy(result_file, os.path.join(Ouput_folder, file_name))
                        break  #We have a hit
    return

# Retrieve the protein sequence and the localisation
def Retrieve_sequence(assembly, sequence_id):
    protein_faa_path = os.path.join(assembly, "genomic.gbff")
    retrieve = []
    for record in SeqIO.parse(protein_faa_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                if sequence_id in feature.qualifiers["protein_id"]:
                    retrieve.append(feature.location)
                    if feature.qualifiers.get(["Translation not available"][0]):
                        print("Translation not available")
                        exit()
                    retrieve.append(feature.qualifiers.get("translation", ["Translation not available"])[0])
                    
                    return retrieve

# Process the file and extract relevant information from the hits.
def process_file(args):

    file_name, Hits_profile, Evalue, Proteome_path,Model_name,Filtered_Hits = args
    # Create a local dictionary to aggregate results
    local_dico = {"bitscore":[], "evalue":[], "id":[],"sequence":[],"assembly":[],"genomic_location":[],"description":[]}  
    result_file = os.path.join(Hits_profile, file_name)
    assembly = file_name.replace("_Profile.hmmer", "")
    assembly = assembly.replace(f"{Model_name}_","")

    assembly_path = os.path.join(Proteome_path, assembly)  
    with open(result_file, 'r') as f:
        
        for queryresult in SearchIO.parse(f, "hmmscan3-domtab"):
            for hit in queryresult.hits:
                # If we have a hit below the evalue
                if hit.evalue <= Evalue:
                    # Copy the result file in the Filtered folder
                    os.makedirs(Filtered_Hits, exist_ok=True)
                    shutil.copy(result_file, os.path.join(Filtered_Hits, file_name))
                    #Get the sequence and the location
                    retrieve = Retrieve_sequence(assembly_path, hit.id)
                    start = 0
                    end = 0
                    for hsp in hit.hsps:
                        start = hsp.query_start
                        end = hsp.query_end
                    for attrib in attributs_hits:
                        if attrib == "sequence":
                            seq = retrieve[1]
                            local_dico[attrib].append(seq[start:end])
                        elif attrib == "assembly":
                            local_dico[attrib].append(assembly)
                        elif attrib == "genomic_location":
                            local_dico[attrib].append(retrieve[0])
                        else:
                            local_dico[attrib].append(getattr(hit, attrib))
    return local_dico

#Filter the results based on evalue using multiprocessing
def Evalue_Filter(Hits_profile, Evalue, Proteome_path,Model_name,Filtered_Hits):
    # get the number of CPUs available
    cpu_number = cpu_count()
    # Create a pool of worker processes
    pool = Pool(processes=cpu_number)  
    # Create a list of arguments to pass to the process function
    args_list = [(file_name, Hits_profile, Evalue, Proteome_path,Model_name,Filtered_Hits) for root, dirs, files in os.walk(Hits_profile) for file_name in files]
    # Use the pool of processes to execute the process function
    results = pool.map(process_file, args_list)
    pool.close()
    pool.join()
    
    # Initialize a global dictionary to aggregate results
    global_dico = {"bitscore":[], "evalue":[], "id":[],"sequence":[],"assembly":[],"genomic_location":[],"description":[]}  
    # Aggregate results from all files
    for result in results:
        for key in global_dico.keys():
            global_dico[key].extend(result[key])
    return global_dico

#Main function to apply the HMMER search in each proteome§
def profiles(Proteome_path,Model_file,Bank,Evalue,Results_path,Model_name):
   

    Hmmer_output_path = f"{Results_path}/All_Profiles/{Bank}"
    Hmmer_hits_output_path = f"{Results_path}/Hits_Profiles/{Bank}"
    #Apply the HMMER search in each proteome
    Hmmer_Search(Model_file, Proteome_path,Hmmer_output_path,Bank,Model_name)
    #Filter the HMMER hits
    Hmmer_Hits_Filter(Hmmer_output_path, Hmmer_hits_output_path)

    Filtered_Hits = f"{Results_path}/Hits_Profiles_Filtered/{Bank}"
    aggregated_results = Evalue_Filter(Hmmer_hits_output_path, Evalue, Proteome_path,Model_name,Filtered_Hits)
    pd_dict = pd.DataFrame.from_dict(aggregated_results)


    pd_dict.to_csv(f"{Results_path}/{Model_name}_Filtered_Sequence_{Bank}.csv",index=False)

    return f"{Results_path}/{Model_name}_Filtered_Sequence_{Bank}.csv"


