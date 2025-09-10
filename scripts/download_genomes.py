import csv
import subprocess
import os

csv_file = "taxid_to_accnum.csv"
output_dir = "downloads"

os.makedirs(output_dir, exist_ok=True)

with open(csv_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        taxid = row["taxid"].strip()
        accession = row["accession_number"].strip()

        if accession == "NOTFOUND" or accession == "":
            print(f"Skipping {taxid} (no accession)")
            continue

        zip_file = os.path.join(output_dir, f"{taxid}.zip")
        folder_name = os.path.join(output_dir, taxid)

        if os.path.exists(folder_name):
            print(f"Already downloaded {taxid}, skipping.")
            continue

        print(f"Downloading {accession} for taxid {taxid}...")

        # Run datasets download
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "gff3,rna,cds,protein,genome,seq-report",
            "--filename", zip_file
        ]

        try:
            subprocess.run(cmd, check=True)

            # unzip
            subprocess.run(["unzip", "-q", "-o", zip_file, "-d", output_dir], check=True)

            # rename ncbi_dataset folder to taxid
            ncbi_folder = os.path.join(output_dir, "ncbi_dataset")
            if os.path.exists(ncbi_folder):
                os.rename(ncbi_folder, folder_name)

            print(f"✓ Saved {taxid} in {folder_name}")

        except subprocess.CalledProcessError as e:
            print(f"Error downloading {accession} for taxid {taxid}: {e}")
