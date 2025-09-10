# Assessing read-stealing (process log)

We want to ensure that reads are not being "stolen" by the contaminant genomes.
Recall that an early step in the processing pipeline is to align raw reads to contaminant genomes (CGs).
These we'll term *rejected reads*, as they are filtered out of the remainder of the pipeline, and not assembled into contigs.
Note that each contaminant genome is indexed by a specific taxid.

JGI provides data relating to this filtering step for each of the ~250 samples.
(Note: For convenience, we'll refer to samples via their unique sample ID, or `sampid` in code.)
In particular, for every sample there is a base coverage file `<sampid>.fastq.gz/<sampid>.sam.pileup.basecov.gz` that provides the number of rejected reads mapped to each nucleotide of every contaminant genome.
Note that every rejected read should be mapped to a unique contaminant genome, as suggested by the fact that the CGs are concatenated into a single file.

An immediate question to consider is whether a suspicious number of rejected reads have been mapped to regions of the contaminant genomes that correspond to our proteins of interest: nar, nir, etc..
To this end, we can assess the depth data contained in the base coverage file for each sample.
Since we are working primarily at the level of KOs, we will use this level of coarse graining.
We need to identify the regions of each genome that are associated to our KOs of interest, then compute the average depth of reads aligned to these regions.

Our approach is outlined as follows:

0. Obtain contaminant genomes and annotations to match the JGI pipeline.
1. Obtain for every KO of interest, a set of representative amino acid sequences (across species) associated to the KO.
2. Use these representative sequences and `diamond` to annotate the genes of each contaminant genome with a KO.
3. For every KO, and every contaminant genome, identify the regions of interest along the genome associated to the KO.
4. For every KO, and every contaminant genome, locate those identified regions of interest on the genome.
5. For every sample, use the sample's coverage file to compute the average depth of reads along each of the identified regions of interest, and assemble this data into tables for downstream analysis.

## Environment setup

```bash
conda create -n bioenv
```

```bash
conda create -n ncbienv
```

## Full pipeline commands

```bash
#~~~  Step 0  ~~~#
conda activate ncbienv
python scripts/download_genomes.py
sh scripts/move_downloaded_data.sh

sh scripts/extract_genome_cds_and_proteins.sh

#~~~  Step 1  ~~~#
python scripts/fetch_ko_info.py -i data/ko_list.txt -o data/ko_info

for ko in $(cat data/ko_list.txt); do 
    echo $ko 
    python scripts/fetch_ko_seqset.py -k ${ko} -o data/kos/${ko} -b 10 --pbar
done

for ko in $(cat data/ko_list.txt); do 
    if [[ -d data/kos/$ko ]]; then 
        echo $ko 
        ./scripts/cluster_ko_seqset.sh $ko
    fi
done

> data/all_kos_rep.faa
for ko in $(ls data/kos); do 
    echo $ko
    cat data/kos/$ko/${ko}_rep.faa >> data/all_kos_rep.faa
done
awk -F' ' '/^>/ {print substr($1,2), $2}' data/all_kos_rep.faa > data/geneid_to_ko.txt

#~~~  Step 2  ~~~#
./scripts/build_diamond_db.sh

dbdir=out/diamond_db
for ko in $(cat data/ko_list.txt); do 
    if [[ -d data/kos/$ko ]]; then 
        echo $ko
        for db in ${dbdir}/*; do
            taxid=$(basename $db .dmnd)
            ./scripts/run_diamond.sh $taxid $ko
        done
    fi
done

./scripts/merge_diamond_results.sh

python scripts/filter_best_matches.py

#~~~  Step 3  ~~~#
./scripts/filter_unique_diamond_hits.sh

#~~~  Step 4  ~~~#
./scripts/run_all_locate_regions_in_genome.sh

#~~~  Step 5  ~~~#
python scripts/compute_region_counts.py

```

## Processing description

### Step 0: Aquiring contaminant genomes and annotations

Download genomes and annotation files using the `download_genome.py` script, which reads from the file `data/taxid_to_accnum.csv`. Then move the downloaded files using the shell script `move_downloaded_data.sh`.

```bash
conda activate ncbienv
python scripts/download_genomes.py
sh scripts/move_downloaded_data.sh
```

Now let's extract from each of the `cds_from_genomic.fna` files a list of the CDS names, and from each of the `protein.faa` files a list of the protein names.

```bash
sh scripts/extract_genome_cds_and_proteins.sh
```

<!-- With these, we can now check to see which genes, based on search keys, are found in each genome.
Running the script `search_genome.sh` will produce a number of csv files, with the results of searching each individual genome for the specified genes.
Once we run this on each individual genome, we want to merge the results into a single table, where the rows correspond to each scaffold, and columns correspond to searches. -->


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

Finally, we can merge all representative sequences into a single file for convenience, and map the sequence identifiers with their assigned KO.

```bash
> data/all_kos_rep.faa
for ko in $(ls data/kos); do 
    echo $ko
    cat data/kos/$ko/${ko}_rep.faa >> data/all_kos_rep.faa
done
awk -F' ' '/^>/ {print substr($1,2), $2}' data/all_kos_rep.faa > data/geneid_to_ko.txt
```

### Step 2: KO annotation of contaminant genomes with `diamond`

We can now search the contaminant genomes for these KO-associated sequences.
We use `diamond` for this.
First, we construct a diamond database per genome. (Recall, we use the taxid to index each genome.)

```bash
./scripts/build_diamond_db.sh
```

This will create and populate a directory `out/diamond_db` with a diamond database file per taxid: `<taxid>.dmnd`.
Once the databases are constructed, we can query the genomes for the representative sequences, as follows:

```bash
taxid=<TAXID>
ko=<KO>
./scripts/run_diamond.sh ${taxid} ${ko}

# For all KOs and all genomes:
dbdir=out/diamond_db
for ko in $(cat data/ko_list.txt); do 
    if [[ -d data/kos/$ko ]]; then 
        echo $ko
        for db in ${dbdir}/*; do
            taxid=$(basename $db .dmnd)
            ./scripts/run_diamond.sh $taxid $ko
        done
    fi
done
```

The resulting files are stored as `out/diamond_res/<ko>/<taxid>_hits.tsv` and have the following structure:

```txt
# Example file: out/diamond_res/K00370/305_hits.tsv

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
The following script constructs this combined file and saves it as `out/dmnd_combined.tsv`.

```bash
./scripts/merge_diamond_results.sh
```

The script `scripts/filter_best_matches.py` takes this combined file, and processes the non-unique maps from region IDs to KOs.
That is, our `diamond` query might associate a particular region of a contaminant genome to multiple KOs, if two different sequences from the KEGG database, associated with different KOs, both sufficiently match to the same contaminant genome region.
In this case, we select the best hit by first minimizing the `evalue` field and then maximizing the `pident` field.
The resulting file is stored as `out/dmnd_combined_best_match.tsv`.

```bash
python scripts/filter_best_matches.py
```

### Step 3: Identifying KO-associated regions within contaminant genomes

From the previous step, we'll construct a list of the unique IDs of sequence hits identified per KO per contaminant genome, and then locate these regions in each of the contaminant genomes.

```bash
# Run via script: ./scripts/filter_unique_diamond_hits.sh

outdir=out/diamond_res_uniq
mkdir -p $outdir
for ko in $(ls out/diamond_res); do
    echo $ko
    for resfile in out/diamond_res/$ko/*.tsv; do
        base=$(basename $resfile .tsv)
        taxid=${base/_hits/}
        echo "Finding unique sequence hits for $ko (taxid=$taxid)"
        cut -f2 "$resfile" | sort | uniq > "${outdir}/${taxid}_res_unique.txt"
    done
done
```

This process is wrapped in the following script:

```bash
./scripts/filter_unique_diamond_hits.sh
```

The resulting output looks as follows:

```txt
# Example file: out/diamond_res_uniq/K00370/305_res_unique.txt

WP_075463358.1
WP_075465178.1
WP_075468174.1
```

Now we have for each KO and each contiminant genome a list of the unique region IDs within the contaminant genome that are similar to at least one of the KO's representative amino acid sequences.

### Step 4: Locating regions of interest within the contaminant genome

We next need to look at the annotation files for our contaminant genomes, find the regions of interest, and identify their start and stop positions within the genome.

The following python script locates the regions along a specified contaminant genome associated to each of a number of KOs.

<!-- ```bash
# Run all with script: ./scripts/run_all_locate_regions.sh

python scripts/locate_regions.py --t <taxid> -k <KO1> <KO2> ... <KO> -a data/genomes/<taxid>/ncbi_data/{accref}/genomic.gff -d out/diamond_res_uniq -o out/identified_regions
``` -->

```bash
# Run all with script: ./scripts/run_all_locate_regions_in_genome.sh

python scripts/locate_regions_in_genome.py -t <taxid> \
    -a data/genomes/<taxid>/ncbi_data/<accref>/genomic.gff \
    -d out/dmnd_combined_best_match.tsv \
    -o out/identified_regions_best_match
```

Running the script `scripts/run_all_locate_regions_in_genome.sh` will result in a number of files in the directory `out/identified_regions`, each corresponding to a contaminant genome.
As an example:

```txt
# Example file: out/identified_regions/identified_regions_305.tsv

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

Our goal now is to compute for each sample the average depth of rejected reads mapped to every identified region of interest of the contaminant genomes.
The following script performs this computation per sample, and merges the results for every contaminant genome into a single file.

```bash
python scripts/compute_region_counts.py
```

The resulting file is of the form shown below.
Each row corresponds to a region of interest identified within one of the contaminant genomes, and associated to a KO.
It provides the start and stop position on the contaminant genome, the sequence length, average read depth, and max depth.

```text
# Example file: out/regions/coverage_Soil3_CE_239_-7_CHL_T9.csv

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
# Example lines from file out/regions/coverage_Soil3_CE_239_-7_CHL_T9.csv

267608,K15864,WP_011004281.1,1244384,1244995,611,1.0,1.0,ID=cds-WP_011004281.1...
267608,K15864,WP_043876895.1,1884433,1885926,1493,10.969859343603483,30.0,"ID=...
267608,K00372,WP_011001989.1,2228951,2231302,2351,10.161207996597193,56.0,"ID=...
```

For example, the lines indicate that for this specific sample, taxon 267608 had an average of 10.97 reads mapped to the 1493-length gene associated with the ID `WP_043876895.1`, which was associated to the KO K15864.

The question now is how one should go about compiling this data.
We want to show samples individually, and elucidate how the rejected reads map to the contaminant genomes.

