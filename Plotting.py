from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import subprocess
import os
import csv
import matplotlib.colors
cluster_to_representative = {} # Mapping between assembly and their representative
representatives = [] # List of representatives
protein_counts_per_rep = defaultdict(lambda: defaultdict(list)) # Store the evalue for each representative and protein

#Convert the DataFrame to FASTA format
def df_to_fasta(df, fasta_file_path):
    # Filter DataFrame to include only assemblies that appear once
    unique_assembly = df['assembly'].value_counts()
    unique_occurrence_assemblies = unique_assembly[unique_assembly == 1].index
    filtered_df = df[df['assembly'].isin(unique_occurrence_assemblies)]
    
    # Write filtered DataFrame to FASTA file
    with open(fasta_file_path, "w") as fasta_file:
        for _, row in filtered_df.iterrows():
            # Writing the FASTA format: ">id\nsequence\n"
            fasta_file.write(f">{row['assembly']}\n{row['sequence']}\n")

# Sort the clusters because the representative is not always in the first line
def sort_clusters(Clusters_path):
    if not os.path.exists(Clusters_path):
        print(f"File not found: {Clusters_path}")
        exit()
    with open(Clusters_path, 'r') as file:
        content = file.read()
    sorted_clusters = []
    clusters = content.split('>Cluster')
    for cluster in clusters[1:]:  
        lines = cluster.strip().split('\n')
        cluster_header = f">Cluster{lines[0]}"
        members = lines[1:]
        # Find the member ending with '*' and place it first after the cluster header
        starred_member = None
        for member in members:
            if member.strip().endswith('*'):
                starred_member = member
                members.remove(member)
                break
        sorted_cluster = [cluster_header, starred_member] + members
        sorted_clusters.append('\n'.join(sorted_cluster))
    sorted_content = '\n'.join(sorted_clusters)
    with open(Clusters_path, 'w') as output_file:
        output_file.write(sorted_content)
     

#Cd-hit to find representative
def Clustering(FilteredSeq,Model_name,Bank,Results_Path):
    df = pd.read_csv(FilteredSeq)
    Fasta_file_path = f"{Results_Path}/{Model_name}_{Bank}_Fasta_file.fasta"
    df_to_fasta(df,Fasta_file_path)

    cdhit_command = ["cd-hit", "-i", f"{Fasta_file_path}","-M","unlimited", "-o", f"{Results_Path}/{Model_name}_{Bank}_clusters", "-c", "0.90", "-n", "5"]
    subprocess.run(cdhit_command, check=True)
    sort_clusters(f"{Results_Path}/{Model_name}_{Bank}_clusters.clstr")
    return f"{Results_Path}/{Model_name}_{Bank}_clusters"

# Get the number of clusters
def Cluster_Number(Clustered_file):
    with open(Clustered_file, "r") as file:
        for line_number, line in enumerate(reversed(file.readlines()), start=1):
            if "Cluster" in line:
                return int(line.split("Cluster")[1])

    return 0

# Add all the non-unique assemblies to the end of the file
def Complete_Clstr(Clustered_file,FilteredSeq):
    df = pd.read_csv(FilteredSeq)

    # Filter DataFrame to include only non-unique assemblies
    unique_assemblies = df['assembly'].value_counts()
    non_unique_assemblies = unique_assemblies[unique_assemblies > 1].index
    non_unique_df = df[df['assembly'].isin(non_unique_assemblies)]
    
    # Compute the length of the sequences
    #seq_lengths = [len(seq) for seq in non_unique_df['sequence']]
    written_assemblies = set()
    # Append non-unique assemblies to the end of the file
    with open(Clustered_file, "a") as cluster_file:
        if (os.stat(Clustered_file).st_size != 0) :
            cluster_file.write("\n")  # Add a new line before appending non-unique assemblies
        current_cluster = Cluster_Number(Clustered_file)  # Current cluster number
        for assembly in non_unique_df['assembly']:
            if assembly not in written_assemblies:
                current_cluster += 1
                cluster_file.write(f">Cluster{current_cluster}\n")
                cluster_file.write(f"0\t{0}aa, >{assembly}... *\n")
                written_assemblies.add(assembly)



# Parsing clusters file to get representatives and their clusters members
def Parsing(Clusters_path):
    global representatives 
    global cluster_to_representative

    current_reprentative = ""
    if (os.stat(f"{Clusters_path}.clstr").st_size == 0):
        return
    with open(f"{Clusters_path}.clstr", 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                True
            else:
                member = line.split('>')[1].split('...')[0]
                if '*' in line:
                    representatives.append(member)
                    current_reprentative = member
                cluster_to_representative[member] = current_reprentative


    return



# Function to create a matrix with representatives and mean e-values for each unique protein
def Representative_Matrix(Filtered_Sequence):
    global protein_counts_per_rep
    global representatives
    unique_proteins = set()

    with open(Filtered_Sequence, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            try:
                _, evalue, protein_id, _, assembly,_,_ = row
                evalue = float(evalue)  
                unique_proteins.add(protein_id)  # Add to our set of unique proteins
                if assembly in cluster_to_representative:
                    rep = cluster_to_representative[assembly]
                    # Append evalue to the list for this protein_id and representative
                    protein_counts_per_rep[rep][protein_id].append(evalue)
            except ValueError:
                # Catch error for seqeuence in mutiple parts
                print(f"Skipping malformed line: {row}")
                continue

    # matrix with representatives and mean e-values for each unique protein
    matrix = [["Representative"] + list(unique_proteins)]
    # Populate the matrix with mean e-values
    for rep in representatives:
        row = [rep]
        for protein_id in unique_proteins:
            if protein_id in protein_counts_per_rep[rep]:
                #mean_evalue = np.mean(protein_counts_per_rep[rep][protein_id])
                
                mean_evalue = np.mean(protein_counts_per_rep[rep][protein_id])
                
                row.append(mean_evalue)
            else:
                row.append(-0.1)  # Use -0.1 if there are no e-values for this protein_id and rep
        matrix.append(row)
    
    return matrix

# Convert matrix into DataFrame 
def Plot_heatmap(matrix,Path_Heatmap,Model_name,Bank):
    df_matrix = pd.DataFrame(matrix[1:], columns=matrix[0])
    # To have the Reprensentative on the x axis
    df_matrix_transposed = df_matrix.set_index('Representative')
    # if no hits --> no plot
    if len(df_matrix_transposed.index) <=0 or len(df_matrix_transposed.columns) <= 0:
        print("No hits --> No HeatMap generated")
        return 0
    # if one hit --> heatmap
    elif len(df_matrix_transposed.index) == 1 or len(df_matrix_transposed.columns) == 1:
        num_rows, num_columns = df_matrix_transposed.shape
        figsize_x = num_columns * 0.25
        figsize_y = num_rows * 0.25

        # Ensure the figsize is within practical limits 
        figsize_x = max(figsize_x, 20) # Max width  between 20 inches and dynamic size
        figsize_y = max(figsize_y, 20) # Max height between 20 inches and dynamic size
        sns.set_theme(font_scale=2)
        plt.figure(figsize=(figsize_x, figsize_y))
        sns.heatmap(df_matrix_transposed, annot=True, cmap="coolwarm",square=True,linecolor="white",linewidths = 0.5,cbar = False,annot_kws={"size": 10})
        plt.title(f"Heatmap for {Model_name}_{Bank}")
        plt.ylabel("Representative")
        plt.xlabel("Protein ID")
        plt.tight_layout()

        plt.savefig(Path_Heatmap)
        plt.close()
        return 1
    # if more than one hit --> clustermap
    else:
        #Compute the figsize of the figure so all the data is visible
        #get the min evalue from the matrix
        v_max = df_matrix_transposed.max().max()
        v_min = df_matrix_transposed.min().min()
        num_rows, num_columns = df_matrix_transposed.shape
        figsize_x = num_columns * 0.25
        figsize_y = num_rows * 0.25
        # Ensure the figsize is within practical limits 
        figsize_x = max(figsize_x, 20) # Max width or 20 
        figsize_y = max(figsize_y, 10) # Max height or 20 

        if(num_rows < 35 or num_columns <35 ):
            sns.set_theme(font_scale=2)
        else:
            sns.set_theme(font_scale=1.3)
        ax = sns.clustermap(df_matrix_transposed,vmax = v_max, vmin = v_min, cmap = "coolwarm", annot=False, figsize=(figsize_x, figsize_y),xticklabels=True, yticklabels=True, square=True,robust = True,linecolor="white",linewidths = 0.5,cbar_kws={
                     "label": "Value scale",
                     "ticks": [v_min, v_max],
                     "format": "%.0f"  # Formatting ticks to show as integers
                 })
        cax = ax.ax_heatmap.collections[0].colorbar
        cax.set_ticks([v_min, v_max])
        cax.ax.set_yticklabels([v_min, v_max])
        ax.figure.suptitle(f"Heatmap for {Model_name}_{Bank}") 


        
        
        plt.savefig(Path_Heatmap)
        plt.close()
        return 2



# Function to extract the genomic location from the filtered sequences csv
def Genomic_Location(Filtered_sequences):
    # Extract the numeric parts using regular expressions
    numbers = Filtered_sequences.split('[')[1].split(']')[0].replace('<', '').replace('>', '').split(':')
    # Convert to integers
    start, end = map(int, numbers)
    return start, end

def Plot_Location(df, Model_name, Bank,Results_Path):
    Output_folder = f"{Results_Path}"
    os.makedirs(Output_folder, exist_ok=True)
    df = df.sort_values(by='start', ascending=True)


    # Adjust figure height dynamically based on the number of proteins
    unique_protein = df["id"]
    fig_height = max(10, len(unique_protein) * 0.5)  # Ensuring a minimum height of 10

    plt.figure(figsize=(15, fig_height))

    # Mapping the unique assemblies to a y-coordinate
    protein_names = df["id"].unique()
    protein_y_coord_map = {assembly: i for i, assembly in enumerate(protein_names)}
    
    # Plot the genomic intervals as horizontal lines
    for index, row in df.iterrows():
        y = protein_y_coord_map[row["id"]]
        plt.hlines(y, xmin=row['start'], xmax=row['end'], color='darkblue', lw=2)
    
    # Adding markers at the start and end points of the intervals
        plt.plot(row['start'], y, marker='o', color='darkblue')
        plt.plot(row['end'], y, marker='o', color='darkblue')

    # Set y-ticks to the center of the intervals
    plt.yticks(ticks=range(len(protein_names)), labels=protein_names)
    plt.tick_params(axis="both",labelsize=15)
    
    
    plt.xlabel('Genomic Position')
    plt.ylabel('Proteins',rotation=90)
    plt.title(f"Genomic Intervals for {Model_name}",fontdict={'fontsize': 16})
    #plt.title(f"Genomic Intervals for {df['id'].iloc[0]}")
    plt.tight_layout()
    plt.savefig(f"{Output_folder}/{Model_name}_{Bank}_Genomic_Location_Plot.png", bbox_inches='tight')

    plt.close()
    return


def Plot(Results_Path,Filtered_sequences,Model_name,Bank):
    global cluster_to_representative, representatives, protein_counts_per_rep

    # Reset global variables otherwise data leak
    cluster_to_representative = {}
    representatives = []
    protein_counts_per_rep = defaultdict(lambda: defaultdict(list))

    Clustered = Clustering(Filtered_sequences,Model_name,Bank,Results_Path)
    Complete_Clstr(f"{Clustered}.clstr",Filtered_sequences)
    Parsing(Clustered)
    matrix = Representative_Matrix(Filtered_sequences)
    #write the matrix as a csv file
    Path_matrix = f"{Results_Path}/{Model_name}_{Bank}_Mean_Evalue_Matrix.csv"
    df_matrix = pd.DataFrame(matrix[1:], columns=matrix[0])
    df_matrix.to_csv(Path_matrix, index=False)

    Path_Heatmap = f"{Results_Path}/{Model_name}_{Bank}_HeatMap.png"
    Number_hits = Plot_heatmap(matrix,Path_Heatmap,Model_name,Bank)
    if Number_hits >=1:
        df = pd.read_csv(Filtered_sequences)    
    # Applying the Genomic_Location function and creating new columns for start and end
        df[['start', 'end']] = df['genomic_location'].apply(
            lambda x: pd.Series(Genomic_Location(x))
        )
        
        Plot_Location(df,Model_name,Bank,Results_Path)
    return
