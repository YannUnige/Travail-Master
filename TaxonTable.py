import os
import pandas as pd
from Bio import SeqIO

# Créer une liste pour stocker les données
data = []

def get_organism_details_from_gbff(file_path):
    details = []
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, "genbank"):
            organism_name = record.annotations.get('organism')
            taxon_id = None

            # Search for the taxon number in the features
            for feature in record.features:
                if feature.type == "source":
                    for xref in feature.qualifiers.get('db_xref', []):
                        if xref.startswith('taxon:'):
                            taxon_id = xref.split(':')[1]
                            break
                    if taxon_id:
                        break

            details.append(organism_name)
            details.append(taxon_id)
            # Optionally, you can return here if you only need the first record
            # return organism_name, taxon_id

            return details  # Returns a list of tuples with (organism_name, taxon_id)

# Fonction pour extraire les informations du fichier GenBank
def extract_gb_info(file_path):
    #record = SeqIO.read(file_path, "genbank")
    organism_name = ""
    taxonomic_number = ""
    organism_name = get_organism_details_from_gbff(file_path)

    # Extraction des données depuis les features
    """ for feature in record.features:
        if feature.type == "source":
            if 'organism' in feature.qualifiers:
                organism_name = feature.qualifiers['organism'][0]
            if 'db_xref' in feature.qualifiers:
                tax_refs = feature.qualifiers['db_xref']
                for ref in tax_refs:
                    if ref.startswith("taxon:"):
                        taxonomic_number = ref.split(":")[1]
                        break
            break """
    
    return organism_name[0], organism_name[1]

# Fonction pour parcourir les dossiers et extraire les informations
def extract_info(path, code_prefix):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.gbff'):
                file_path = os.path.join(root, file)
                # Extraire le code
                code = root.split('/')[-1]
                if code_prefix:
                    code = code_prefix + code

                # Appel de la fonction pour extraire les informations
                organism_name, taxonomic_number = extract_gb_info(file_path)
                
                # Ajouter à la liste
                data.append([code, organism_name, taxonomic_number])

# Extraire les informations du premier chemin
extract_info('Data/Genome/Megaviricetes/ncbi_dataset/data', None)

# Extraire les informations du deuxième chemin
extract_info('Data/Genome/Megaviricetes_NucleotideData', None)

# Créer un DataFrame
df = pd.DataFrame(data, columns=['Code', 'Nom de l\'organisme', 'Numéro taxonomique'])

# Enregistrer le DataFrame dans un fichier Excel
df.to_excel('genomes_info.xlsx', index=False)
