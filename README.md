# Process Log: Assessing read-stealing

We want to ensure that reads are not being "stolen" by the contaminant genomes.
Recall that an early step in the processing pipeline is to align raw reads to contaminant genomes (CGs).
These we'll term *rejected reads*, as they are filtered out of the remainder of the pipeline, and not assembled into contigs.
Note that each contaminant genome is indexed by a specific taxonomic ID (taxid).

JGI provides data relating to this filtering step for each of the ~250 samples.
(Note: For convenience, we'll refer to samples via their unique sample ID, or `sampid` in code.)
In particular, for every sample there is a base coverage file `<sampid>.fastq.gz/<sampid>.sam.pileup.basecov.gz` that provides the number of rejected reads mapped to each nucleotide of every contaminant genome.
Note that every rejected read should be mapped to a unique contaminant genome, as suggested by the fact that the CGs are concatenated into a single file.

An immediate question to consider is whether a suspicious number of rejected reads have been mapped to regions of the contaminant genomes that correspond to our proteins of interest: nar, nir, nap, etc..
To this end, we will assess the depth data contained in the base coverage file for each sample.
Since we are working primarily at the level of KOs, we will use this level of coarse graining.
We need to identify the regions of each genome that are associated to our KOs of interest, then compute the average depth of reads aligned to these regions.

Our approach is outlined as follows:

0. Obtain contaminant genomes and annotations to match the JGI pipeline.
1. Obtain for every KO of interest a set of representative amino acid sequences (across species).
2. Use these representative sequences and `diamond` to find similar genes within each contaminant genome.
3. Associate each region of interest to a unique KO by identifying the best diamond hit.
4. Locate the identified regions of interest within contaminant genomes.
5. Compute average and max read depth at each region of interest.

## Environment setup

```bash
conda env create -n <env-name> -f environment.yml
conda activate <env-name>
```

## Full pipeline commands

The full processing pipeline is contained in and can be run by executing the following script:

```bash
conda activate <env-name>
./scripts/run_full_pipeline.sh <outdir> <pident_thresh> <evalue_thresh>
```

## Processing steps

### Step 0: Aquiring contaminant genomes and annotations

We first download genomes and annotation files using the `download_genome.py` script, which reads from the file `data/taxid_to_accnum.csv`.
Then, we organize the downloaded files using the shell script `move_downloaded_data.sh`.

```bash
conda activate ncbienv
python scripts/download_genomes.py
sh scripts/move_downloaded_data.sh
```

#### Caveats

##### Caveat 1: E. coli spike-in

The spike-in genomes must be handled slightly differently, before proceeding.
For the first spike-in, a strain of E. coli, we will use the genome and annotation provided by Kiseok.
We store this in the `data/genome` directory, mimicking the structure of the NCBI downloads:

```text
genomes/
└── 1/  <--  arbitrary taxid for custom E. coli strain
    └── ncbi_data/
       └── GCF_000000000.1/  <-- arbitrary name to match structure
           ├── cds_from_genomic.fna  <-- Lee_A8Q_1_polished.ffn
           ├── genomic.gff  <-- Lee_A8Q_1_polished.gff3
           └── protein.faa  <-- Lee_A8Q_1_polished.faa
```


##### Caveat 2: Parabacteroides spike-in

For the second spike-in of Parabacteroides sp. MSK.9.14, we find the taxid 2849180 and acquire the genome annotation from NCBI, as for the other contaminant genomes.
However, in the basecoverage files, this genome is referenced through multiple scaffolds instead of just a single scaffold, as is the case for the other contaminant genomes.
If we check the coverage files for this spike-in, we see something like

```txt
# zcat < 53004.7.525186.AATTAGATTG-GCATTTAGCG.sam.pileup.basecov.gz | grep "NZ_JAHPXO010000???.1_Parabacteroides_sp"

NZ_JAHPXO010000001.1_Parabacteroides_sp._MSK.9.14   0   24
NZ_JAHPXO010000001.1_Parabacteroides_sp._MSK.9.14   1   24
NZ_JAHPXO010000001.1_Parabacteroides_sp._MSK.9.14   2   24
...
NZ_JAHPXO010000001.1_Parabacteroides_sp._MSK.9.14   336656   17
NZ_JAHPXO010000001.1_Parabacteroides_sp._MSK.9.14   336657   16
NZ_JAHPXO010000002.1_Parabacteroides_sp._MSK.9.14   0   13
NZ_JAHPXO010000002.1_Parabacteroides_sp._MSK.9.14   1   13
NZ_JAHPXO010000002.1_Parabacteroides_sp._MSK.9.14   2   13
...
NZ_JAHPXO010000002.1_Parabacteroides_sp._MSK.9.14   258623   10
NZ_JAHPXO010000002.1_Parabacteroides_sp._MSK.9.14   258624   10
```

For the other contaminant genomes, as well as the E. coli spike-in, we store the base coverage information in numpy archives `data/coverage_arrays/coverage_arrays_<sample_id>.npz` where the archive has a key for every contaminant genome scaffold name.
For this spike-in, however, we store the basecoverage information for each individual scaffold (1 through 153) as `data/coverage_arrays_parabacteroides/coverage_arrays_<sample_id>.npz`.
The differences are demonstrated below:

```python
import numpy as np
arrs = np.load("data/coverage_arrays/coverage_arrays_Soil3_CE_239_-7_CHL_T9.npz")
arrs_para = np.load("data/coverage_arrays_parabacteroides/coverage_arrays_Soil3_CE_239_-7_CHL_T9.npz")

print("standard contaminant coverage keys:")
for k in list(arrs.keys())[0:5]:
    print(k)
print("parabacteroides coverage keys:")
for k in list(arrs_para.keys())[0:5]:
    print(k)
```

Running this code chunk will result in

```txt
standard contaminant coverage keys:

Lee_A8Q_1_Ecoli_contig_1_polypolish
Ralstonia_solanacearum_strain_KACC_10722
CP013959.1_Staphylococcus_aureus_strain_V605
NC_003902.1_Xanthomonas_campestris_pv._campestris_str._ATCC_33913
NC_003295.1_Ralstonia_solanacearum_GMI1000
...

parabacteroides coverage keys:

NZ_JAHPXO010000099.1_Parabacteroides_sp._MSK.9.14
NZ_JAHPXO010000100.1_Parabacteroides_sp._MSK.9.14
NZ_JAHPXO010000101.1_Parabacteroides_sp._MSK.9.14
NZ_JAHPXO010000102.1_Parabacteroides_sp._MSK.9.14
NZ_JAHPXO010000103.1_Parabacteroides_sp._MSK.9.14
...
```

So, in order to handle this, we run the script

```bash
./scripts/modify_parabact_gff.sh
```

which will produce a file `data/genomes/2849180/protein_to_scaffold.tsv` that maps the protein names to their scaffold ID, and then modify the `protein.faa` file so that the new protein names are of the form `<gene_id>:<scaffold_id>`.
This will eventually allow us to obtain the specific scaffold key from a protein name.


### Continuing on

Next, we extract from each of the `cds_from_genomic.fna` files a list of the CDS names, and from each of the `protein.faa` files a list of the protein names.

```bash
sh scripts/extract_genome_cds_and_proteins.sh
```

### Step 1: Downloading KO sets

Our KOs of interest are listed in the file `data/ko_list.txt`.
For each KO, we can first obtain some basic information from the KEGG database, using the script `scripts/fetch_ko_info.py`.

```bash
python scripts/fetch_ko_info.py -i data/ko_list.txt -o data/ko_info
```

This will populate the output directory `data/ko_info` with a text file for each KO.

Since KOs may encapsulate proteins that are not necessarily close in sequence space, we will acquire all sequences mapped to our KOs of interest from the KEGG database.
To do this, we apply the following script to each KO of interest.

```bash
ko=<KO>
python scripts/fetch_ko_seqset.py -k ${ko} -o data/kos/${ko} -b 10 --pbar

# For all KOs:
for ko in $(cat data/ko_list.txt); do 
    echo $ko 
    python scripts/fetch_ko_seqset.py -k ${ko} -o data/kos/${ko} -b 10 --pbar
done
```

This will download in a .faa file all protein sequences found in the KEGG database corresponding to the specified KO.
Next, with the following script we cluster these files for each KO using the `cd-hit` tool.

```bash
ko=<KO>
./scripts/cluster_ko_seqset.sh ${ko}

# For all KOs:
for ko in $(cat data/ko_list.txt); do 
    if [[ -d data/kos/$ko ]]; then 
        echo $ko 
        ./scripts/cluster_ko_seqset.sh $ko
    fi
done
```

This script will produce files `<ko>_rep.faa` and `<ko>_rep.faa.clstr` with the directory `data/kos/<ko>`.
The former is a representative subset of the input `<ko>.faa` file.
The latter specifies the clusters of sequences associated to each representative.
Taking stock, we now have for each KO a set of representative amino acid sequences associated with that KO.
We can merge all representative sequences into a single file for convenience, and map the sequence identifiers (`seqid`) to their assigned KO.

```bash
> <outdir>/all_kos_rep.faa
for ko in $(ls data/kos); do 
    echo $ko
    cat data/kos/$ko/${ko}_rep.faa >> <outdir>/all_kos_rep.faa
done
awk -F' ' '/^>/ {print substr($1,2), $2}' <outdir>/all_kos_rep.faa > <outdir>/seqid_to_ko.txt
```

The script `scripts/build_ko_info.py` produces a tsv file `data/ko_information.tsv` with information compiled from the fetched KO files.

```bash
python scripts/build_ko_info_file.py -i data/ko_list.txt -d data/ko_info -o data/ko_information.tsv
```

### Step 2: KO annotation of contaminant genomes with `diamond`

We can now search the contaminant genomes for the KO-associated sequences using `diamond`.
First, we construct a diamond database per genome. (Recall, we use the taxid to index each genome.)

```bash
./scripts/build_diamond_db.sh <outdir>
```

This will create and populate a directory `<outdir>/diamond_db` with a diamond database file per taxid: `<taxid>.dmnd`.
Once the databases are constructed, we can query the genomes for the representative sequences, as follows:

```bash
taxid=<TAXID>
ko=<KO>
outdir=<OUTDIR>
evalue_thresh=1e-5
pident_thresh=50
./scripts/run_diamond.sh ${taxid} ${ko} ${outdir} ${evalue_thresh} ${pident_thresh}

# For all KOs and all genomes:
dbdir=<outdir>/diamond_db
for ko in $(cat data/ko_list.txt); do 
    if [[ -d data/kos/$ko ]]; then 
        echo $ko
        for db in ${dbdir}/*; do
            taxid=$(basename $db .dmnd)
            ./scripts/run_diamond.sh $taxid $ko $evalue_thresh $pident_thresh
        done
    fi
done
```

The resulting files are stored as `<outdir>/diamond_res/<ko>/<taxid>_hits.tsv` and have the following structure:

```txt
# Example file: <outdir>/diamond_res/K00370/305_hits.tsv

ecos:EC958_1735   WP_075468174.1  65.1  1247  408  9   42  1275  1  1233  0.0  1731
sef:UMN798_1652   WP_075468174.1  63.7  750   250  6   1   738   1  740   0.0  1029
eec:EcWSU1_02361  WP_075468174.1  65.0  1248  410  9   16  1250  1  1234  0.0  1738
enf:AKI40_3370    WP_075468174.1  65.4  1248  406  10  30  1265  1  1234  0.0  1763
...
```

with tab-separated columns

```txt
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

The first column `qseqid` is the ID of the query (species and associated gene ID).
The second column `sseqid` is the ID of the subject sequence identified in the searched contaminant genome.
We want to identify these sequences so that we can examine their average read depth.
For ease of downstream processing, we can concatenate these files into a single file, with additional columns `taxid` and `ko`.
The following script constructs this combined file and saves it as `<outdir>/dmnd_combined.tsv`.

```bash
./scripts/merge_diamond_results.sh
```

Finally, the following script populates a directory `diamond_hits` with the same structure as `diamond_res`, with text files giving the set of seqids found via diamond, for each KO and contaminant genome.

```bash
./scripts/get_diamond_hits.sh <outdir>
```

The resulting output looks as follows:

```txt
# Example file: <outdir>/diamond_hits/K00370/305_res_unique.txt

WP_075463358.1
WP_075465178.1
WP_075468174.1
```

Now we have for each KO and each contiminant genome a list of the unique region IDs within the contaminant genome that are similar to at least one of the KO's representative amino acid sequences.

### Step 3: Associate regions of interest to a unique KO

The script `scripts/filter_top_diamond_hits.py` takes the combined file of `diamond` results, and constructs a mapping from  processes the non-unique maps from seqids to KOs.
That is, our `diamond` query might associate a particular gene in a contaminant genome to multiple KOs, if two different sequences from the KEGG database, associated with different KOs, are both sufficiently close to the same gene.
In this case, we select the best diamond hit by first minimizing the `evalue` field and then maximizing the `pident` field.
The following script performs this filtering step:

```bash
python scripts/filter_top_diamond_hits.py -i <outdir>/dmnd_combined.tsv -o <outdir>/dmnd_combined_top_hits.tsv
```

The resulting file is stored as `<outdir>/dmnd_combined_top_hits.tsv`, and defines a unique mapping from contaminant genes to an assigned KO.


### Step 4: Locating regions of interest within the contaminant genomes

We next need to look at the annotation files for our contaminant genomes, find the regions of interest, and identify their start and stop positions within the genome.
The python script `scripts/locate_regions_in_genome.py` processes one contaminant genome at a time, locating each region listed in the file `dmnd_combined_top_hits.tsv`.

```bash
# Run all with script: ./scripts/run_all_locate_regions_in_genome.sh

python scripts/locate_regions_in_genome.py \
        -t <taxid> \
        -a data/genomes/<taxid>/ncbi_data/<accnum>/genomic.gff \
        -d <outdir>/dmnd_combined_top_hits.tsv \
        -o <outdir>/identified_regions
```

Running the script `scripts/run_all_locate_regions_in_genome.sh` will result in a number of files in the directory `<outdir>/identified_regions`, each corresponding to a contaminant genome.
As an example:

```txt
# Example file: <outdir>/identified_regions/identified_regions_305.tsv

ko      name            start    stop     desc
K15864  YP_005227361.1  3049987  3051327  ID=cds-YP_005227361.1;Dbxref=GenBank:YP_005227361.1...;Name=YP_005227361.1;...
K00372  YP_005224340.1  46224    48638    ID=cds-YP_005224340.1;Dbxref=GenBank:YP_005224340.1...;Name=YP_005224340.1;...
K00372  YP_005224622.1  362708   364387   ID=cds-YP_005224622.1;Dbxref=GenBank:YP_005224622.1...;Name=YP_005224622.1;...
...
K00371  WP_075468172.1  1630550  1632103  ID=cds-WP_075468172.1;Parent=gene-LBM2029_RS22645;...;Name=WP_075468172.1;...
K00370  WP_075463358.1  614765   616861   ID=cds-WP_075463358.1;Parent=gene-LBM2029_RS02845;...;Name=WP_075463358.1;...
K00370  WP_075465178.1  2830618  2832744  ID=cds-WP_075465178.1;Parent=gene-LBM2029_RS12865;...;Name=WP_075465178.1;...
K00370  WP_075468174.1  1632100  1635837  ID=cds-WP_075468174.1;Parent=gene-LBM2029_RS22650;...;Name=WP_075468174.1;...
K04748  WP_075463912.1  1507143  1509458  ID=cds-WP_075463912.1;Parent=gene-LBM2029_RS06925;...;Name=WP_075463912.1;...
...
...
```

The third and fourth columns give the start and stop position of the region.

### Step 5: Computing average read depth

Our goal now is to compute for each sample the average depth of rejected reads mapped to every identified region of interest across the contaminant genomes.
The following script performs this computation per sample, merging results for every contaminant genome into a single file per sample.

```bash
python scripts/compute_region_counts.py \
    -r <outdir>/identified_regions \
    -o <outdir>/ko_expression \
    --coverage_dir data/coverage_arrays
```

The resulting file is of the form shown below.
Each row corresponds to a region of interest identified within one of the contaminant genomes, and associated to a KO.
It provides the start and stop position on the contaminant genome, the sequence length, average read depth, and max depth.

```text
# Example file: <outdir>/ko_expression/coverage_Soil3_CE_239_-7_CHL_T9.csv

taxid,ko,name,start,stop,seq_length,avg_depth,max_depth,desc
267608,K15864,WP_011004281.1,1244384,1244995,611,1.0,1.0,ID=cds-WP_011004281.1...
267608,K15864,WP_043876895.1,1884433,1885926,1493,10.969859343603483,30.0,"ID=...
267608,K00372,WP_011001989.1,2228951,2231302,2351,10.161207996597193,56.0,"ID=...
...
305,K00371,WP_075468172.1,1630550,1632103,1553,0.0,0.0,"ID=cds-WP_075468172.1;...
305,K00370,WP_075463358.1,614765,616861,2096,1.297232824427481,6.0,"ID=cds-WP_...
305,K00370,WP_075465178.1,2830618,2832744,2126,0.21307619943555975,2.0,"ID=cds...
305,K00370,WP_075468174.1,1632100,1635837,3737,0.08081348675408082,1.0,"ID=cds...
305,K04748,WP_075463912.1,1507143,1509458,2315,3.023758099352052,29.0,"ID=cds-...
...
```

We now consider approaches for downstream analysis of these files.


## Downstream analysis

For any given sample, we can now specify a contaminant genome, and query the average depth of reads at any of the regions of interest that were associated with a KO.

```txt
# Example lines from file <outdir>/ko_expression/coverage_Soil3_CE_239_-7_CHL_T9.csv

267608,K15864,WP_011004281.1,1244384,1244995,611,1.0,1.0,ID=cds-WP_011004281.1...
267608,K15864,WP_043876895.1,1884433,1885926,1493,10.969859343603483,30.0,"ID=...
267608,K00372,WP_011001989.1,2228951,2231302,2351,10.161207996597193,56.0,"ID=...
```

For example, the lines indicate that for this specific sample, taxon 267608 had an average of 10.97 reads mapped to the 1493-length gene associated with the ID `WP_043876895.1`, which was associated to the KO K15864.

