import os


def List_Hits_txt(Result_Folder,organisme_name):
    folds_with_hits = {}
    folds_without_hits = {}

    for root, dirs, file in os.walk(Result_Folder):
        model = os.path.basename(root)
        
        if "_" not in model or "CLASS" not in model:  # Skip if naming convention not matched
            continue
        folds_name,class_name = model.split("_CLASS_")
        if "Hits_Profiles_Filtered" in dirs: # If we have hits below evalue for specific model

            if folds_name in folds_with_hits: # if we already have hits for this model
                folds_with_hits[folds_name].append(class_name)
            else:
                folds_with_hits[folds_name] = [class_name] # else add it

        else: # if we don't have hits
            if folds_name in folds_without_hits: 
                folds_without_hits[folds_name].append(class_name)
            else:
                folds_without_hits[folds_name] = [class_name]
                #f.write(model + "\n")
  
    write_to_file(f"{Result_Folder}/{organisme_name}_With_Hits.txt", folds_with_hits)
    write_to_file(f"{Result_Folder}/{organisme_name}_Without_Hits.txt", folds_without_hits)



def write_to_file(filename, data):
    with open(filename, "w") as f:
        for folds_name, classes in data.items():
            f.write(f"{folds_name}\n")
            for class_name in sorted(set(classes)):  # Sort and remove duplicates
                f.write(f"\t- {class_name}\n")
            f.write("\n")


