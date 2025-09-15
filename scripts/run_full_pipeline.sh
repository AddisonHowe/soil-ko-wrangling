#!/usr/bin/env bash

SKIP_GENOME_DOWNLOAD=true
SKIP_KO_DOWNLOAD=true

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 0: Obtain contaminant genomes and annotations          ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 0  ***************"

if [ "$SKIP_GENOME_DOWNLOAD" = true ]; then
    echo "*** Skipping genome downloads!"
else
    echo "*** Downloading genomes..."
    python scripts/download_genomes.py
    sh scripts/move_downloaded_data.sh
fi

echo "*** Extracting genome CDS and proteins..."
sh scripts/extract_genome_cds_and_proteins.sh

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 1: Obtain KO information and representative sequences  ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 1  ***************"

if [ "$SKIP_KO_DOWNLOAD" = true ]; then
    echo "*** Skipping KO downloads!"
else
    echo "*** Fetching KO information..."
    python scripts/fetch_ko_info.py -i data/ko_list.txt -o data/ko_info

    echo "*** Fetching KO representative sequences..."
    for ko in $(cat data/ko_list.txt); do 
        echo $ko 
        python scripts/fetch_ko_seqset.py -k ${ko} -o data/kos/${ko} -b 10 --pbar
    done
fi

echo "*** Clustering KO representative sequences..."
for ko in $(cat data/ko_list.txt); do 
    if [[ -d data/kos/$ko ]]; then 
        echo $ko 
        ./scripts/cluster_ko_seqset.sh $ko
    fi
done

echo "*** Merging all representative sequences across KOs into one file..."
> data/all_kos_rep.faa
for ko in $(ls data/kos); do 
    echo $ko
    cat data/kos/$ko/${ko}_rep.faa >> data/all_kos_rep.faa
done
awk -F' ' '/^>/ {print substr($1,2), $2}' data/all_kos_rep.faa > data/seqid_to_ko.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 2: Run diamond to find near sequences                  ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 2  ***************"

echo "*** Building diamond database..."
./scripts/build_diamond_db.sh

echo "*** Running diamond..."
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

echo "*** Merging diamond results..."
./scripts/merge_diamond_results.sh

echo "*** Building lists of diamond hits..."
./scripts/get_diamond_hits.sh

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 3: Associate regions of interest to a unique KO        ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 3  ***************"

echo "*** Identifying top diamond hits..."
python scripts/filter_top_diamond_hits.py

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 4: Locate regions of interest in contaminant genomes   ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 4  ***************"

echo "*** Locating regions..."
./scripts/run_all_locate_regions_in_genome.sh

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 5: Compute avg read depth at each region of interest   ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 5  ***************"

echo "*** Computing region counts..."
python scripts/compute_region_counts.py \
    -r out/identified_regions \
    -o out/ko_expression \
    --coverage_dir data/coverage_arrays \
    --para_coverage_dir data/coverage_arrays_parabacteroides
