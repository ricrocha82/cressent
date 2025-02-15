import os
import pandas as pd
from Bio import Entrez

# some taxonomies couldn't be found using the query from NCBI
# the accession numbers were acquired from the ICTV reports

family_file='/fs/project/PAS1117/ricardo/ssDNA_tool/build_db/accession_not_found_refseq.csv'
dir_to_download='/fs/scratch/Sullivan_Lab/Ricardo/ssDNA_db/'

# Set your email (Entrez requires an email for API access)
Entrez.email = 'pavan.4@osu.edu'

# Create output directory if it doesn't exist
os.makedirs(dir_to_download, exist_ok=True)

# Read CSV using pandas
df = pd.read_csv(family_file)

# Iterate over each row
for _, row in df.iterrows():
    project, family = row[0], row[1]  # Read Accession and Family

    print(f"Processing: Accession = {project}, Family = {family}")

    # Create a directory for the family
    family_dir = os.path.join(dir_to_download, family)
    os.makedirs(family_dir, exist_ok=True)

    # Output files
    # out_1 = os.path.join(family_dir, f"{project}.txt")  # Stores protein IDs
    out_2 = os.path.join(family_dir, f"{family}.fa")    # Stores protein sequences

    # Search for nucleotide records
    try:
        search_handle = Entrez.esearch(db="nuccore", term=project, retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            print(f"No results found for {project}")
            continue

        nucleotide_id = search_results["IdList"][0]
        print(f"Found nucleotide ID: {nucleotide_id}")

        # Fetch feature table (ft) to extract protein IDs
        fetch_handle = Entrez.efetch(db="nuccore", id=nucleotide_id, rettype="ft", retmode="text")
        feature_table = fetch_handle.read()
        fetch_handle.close()

        # Extract protein IDs from the feature table
        protein_ids = []
        for line in feature_table.splitlines():
            if "protein_id" in line:
                parts = line.split()
                if len(parts) > 1:
                    protein_id = parts[1].strip('"')
                    protein_ids.append(protein_id)

        if not protein_ids:
            print(f"No protein IDs found for {project}")
            continue

        # Retrieve and save protein sequences
        with open(out_2, "a") as fasta_file:
            for protein_id in protein_ids:
                print(f"Downloading: {protein_id}")
                try:
                    fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                    fasta_data = fetch_handle.read()
                    fetch_handle.close()
                    fasta_file.write(fasta_data)
                except Exception as e:
                    print(f"Error fetching protein {protein_id}: {e}")

    except Exception as e:
        print(f"Error processing {project}: {e}")