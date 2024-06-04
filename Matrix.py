import os
import pandas as pd
from openpyxl.utils import get_column_letter

# Generate a matrix that contains occurences of hits for each model and each genome
def Matrix_Hits(genome_dir,model_dir,results_dir,bank,organisme_name):

    genome_names = [name for name in os.listdir(genome_dir) if os.path.isdir(os.path.join(genome_dir, name)) and name.startswith(bank)]

    model_files = os.listdir(model_dir) 
    data_structure = {}

    for model_file in model_files:

        # Extract fold name and class name

        fold_name, class_name = model_file.split("_CLASS_")
        class_name = class_name.replace(".hmm", "")

        # Add the fold name and class name to the data structure
        if fold_name not in data_structure:
            data_structure[fold_name] = {}

        # Add the genome name to the data structure and set its value to 0
        data_structure[fold_name][class_name] = {genome_name: 0 for genome_name in genome_names}
    
    for fold_name in data_structure:

        for class_name in data_structure[fold_name]:
            # Recreate the hits_filtered directory for each model
            dir_path = os.path.join(f"{results_dir}",f"{fold_name}_CLASS_{class_name}")
            hits_dir = os.path.join(dir_path, f"Hits_Profiles_Filtered/{bank}")
            if os.path.exists(hits_dir):
                for genome_name in genome_names:
                    expected_file = f"{fold_name}_CLASS_{class_name}_{genome_name}_Profile.hmmer"
                    if os.path.exists(os.path.join(hits_dir, expected_file)):
                        data_structure[fold_name][class_name][genome_name] = 1

    # Convert data_structure to a pandas DataFrame
    for fold_name in list(data_structure.keys()):
        for class_name in list(data_structure[fold_name].keys()):
            # If all values for this class are 0, delete the class_name key
            if all(value == 0 for value in data_structure[fold_name][class_name].values()):
                del data_structure[fold_name][class_name]
    df = pd.DataFrame.from_dict({(i,j): data_structure[i][j] 
                                for i in data_structure.keys() 
                                for j in data_structure[i].keys()},
                                orient="columns")

    # Save DataFrame to Excel
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    with pd.ExcelWriter(f"{results_dir}/{organisme_name}_{bank}_Hits_matrix.xlsx", engine="xlsxwriter") as writer:
        df.to_excel(writer, sheet_name="Sheet1")
        workbook  = writer.book
        worksheet = writer.sheets['Sheet1']
        format1 = workbook.add_format({'bg_color': '#C6EFCE'})

        # Determine the last column letter and row number otherwise conditional formatting won't apply to the whole sheet

        (max_row, max_col) = df.shape
        col_letter = get_column_letter(max_col + 1)  # Adjust for index column
        max_row += 3  # Adjust for header row

        # Construct the dynamic cell range for conditional formatting
        cell_range = f"B2:{col_letter}{max_row}"

        worksheet.conditional_format(cell_range, {"type": "cell","criteria": "=","value": 1,"format": format1}) 
    
    
