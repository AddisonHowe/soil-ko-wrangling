#!/usr/bin/env bash

SKIP_GENOME_DOWNLOAD=true
SKIP_KO_DOWNLOAD=true

EVALUE_THRESH_DEFAULT=1e-5

if [ "$#" -eq 3 ]; then
    outdir=$1
    pident_thresh=$2
    evalue_thresh=$3
elif [ "$#" -eq 2 ]; then
    outdir=$1
    pident_thresh=$2
    evalue_thresh=$EVALUE_THRESH_DEFAULT
else
    echo "Usage: $0 outdir pident_thresh [evalue_thresh=1e-5]"
    exit 1
fi

# Check if pident_thresh is a file
if [[ -e $pident_thresh ]]; then
    USE_PIDENT_FPATH=true
    pident_fpath=$pident_thresh
    pident_thresh=""
    echo "*** Reading pident values from file $pident_fpath"
    # Enable associative arrays
    declare -A pident_table

    # Read the file and populate the associative array
    while read -r col1 col2 col3; do
        key="${col1}_${col2}"   # combine first two columns into a key
        pident_table["$key"]="$col3"
    done < $pident_fpath
else
    USE_PIDENT_FPATH=false
    pident_fpath=""
fi

mkdir -p $outdir

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
> ${outdir}/all_kos_rep.faa
for ko in $(ls data/kos); do 
    echo $ko
    cat data/kos/$ko/${ko}_rep.faa >> ${outdir}/all_kos_rep.faa
done
awk -F' ' '/^>/ {print substr($1,2), $2}' ${outdir}/all_kos_rep.faa > ${outdir}/seqid_to_ko.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 2: Run diamond to find near sequences                  ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 2  ***************"

echo "*** Building diamond database..."
./scripts/build_diamond_db.sh $outdir

echo "*** Running diamond..."
dbdir=${outdir}/diamond_db
for ko in $(cat data/ko_list.txt); do 
    if [[ -d data/kos/$ko ]]; then 
        printf "%s " "$ko"
        for db in ${dbdir}/*; do
            taxid=$(basename $db .dmnd)
            if [[ $USE_PIDENT_FPATH -eq "true" ]]; then
                lookup_key="${taxid}_${ko}"
                pident_thresh=${pident_table[$lookup_key]}                
            fi
            printf "%s (pi=%.1f) " "$taxid" "$pident_thresh"
            ./scripts/run_diamond.sh $taxid $ko $outdir $evalue_thresh $pident_thresh
        done
        printf "\n"
    fi
done

echo "*** Merging diamond results..."
./scripts/merge_diamond_results.sh $outdir

echo "*** Building lists of diamond hits..."
./scripts/get_diamond_hits.sh $outdir

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 3: Associate regions of interest to a unique KO        ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 3  ***************"

echo "*** Identifying top diamond hits..."
python scripts/filter_top_diamond_hits.py \
    -i ${outdir}/dmnd_combined.tsv \
    -o ${outdir}/dmnd_combined_top_hits.tsv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 4: Locate regions of interest in contaminant genomes   ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 4  ***************"

echo "*** Locating regions..."
./scripts/run_all_locate_regions_in_genome.sh $outdir

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~  Step 5: Compute avg read depth at each region of interest   ~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "**************  Step 5  ***************"

echo "*** Computing region counts..."
python scripts/compute_region_counts.py \
    -r ${outdir}/identified_regions \
    -o ${outdir}/ko_expression \
    --coverage_dir data/coverage_arrays \
    --para_coverage_dir data/coverage_arrays_parabacteroides
