#!/bin/bash

# Base directories
bam_dir="/home/alatypova/At_paper_2025/ont_methylation"
ref_dir="/home/alatypova/At_paper_2025"
modkit_path="/home/alatypova/dist_modkit_v0.4.4_251055f/modkit"

# Loop through barcodes 19 to 31
for barcode in {24..25}; do
# for barcode in 21; do
    # Format barcode with leading zero if necessary
    barcode=$(printf "%02d" $barcode)

    echo "Processing barcode $barcode..."

    # Check if the no_sup file already exists
    no_sup_file="$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.no_sup.bam"
    sorted_bam_file="$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.bam"
    # if [[ -f "$no_sup_file" ]]; then
    #     echo "File $no_sup_file already exists. Skipping Step 0 for barcode $barcode."
    # else
    #     # Step 0: filter out not primary alignments
    #     echo "Running samtools view not primary for barcode $barcode..."
    #     samtools view -F 256 -bo "$no_sup_file" \
    #         "$bam_dir/sorted_SQK-NBD114-96_barcode$barcode.meth.pseudo.bam"
    # fi
    # Step 0: filter out supplementary alignments
    echo "Running samtools view for no prim barcode $barcode..."
    samtools sort -o "$sorted_bam_file" -@ 100 "$bam_dir/SQK-NBD114-96_barcode${barcode}.meth.tair10.bam"

    samtools view -F 2308 -bo "$no_sup_file" \
        "$sorted_bam_file"

    # Step 0.1: get index
    samtools index -@ 100 "$no_sup_file"


    # Step 1: Run modkit pileup
    echo "Running modkit pileup for barcode $barcode... file $no_sup_file"
    $modkit_path pileup --combine-mods -t 100 --filter-threshold C:0.75 \
        "$no_sup_file" \
        "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.pileup.bed"

    # Step 2: Filter the pileup file by 12th column (Nmod >= 3)
    echo "Filtering pileup file for barcode $barcode..."
    awk '{if ($10 >= 3) print $0}' "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.pileup.bed" > \
        "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.pileup.filt.bed"

    # Step 3: Generate motif BED files for CG, CHG, and CHH
    echo "Generating motif BED files for barcode $barcode..."
    $modkit_path motif bed "$ref_dir/TAIR10_1.fna" CG 0 > \
        "$ref_dir/TAIR10_1.CG.bed"
    $modkit_path motif bed "$ref_dir/TAIR10_1.fna" CHG 0 > \
        "$ref_dir/TAIR10_1.CHG.bed"
    $modkit_path motif bed "$ref_dir/TAIR10_1.fna" CHH 0 > \
        "$ref_dir/TAIR10_1.CHH.bed"

    # Step 4: Intersect filtered pileup with motif BED files
    echo "Intersecting filtered pileup with motif BED files for barcode $barcode..."
    bedtools intersect -sorted -a "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.pileup.filt.bed" \
        -b "$ref_dir/TAIR10_1.CG.bed" > \
        "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.CG.bed"
    bedtools intersect -sorted -a "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.pileup.filt.bed" \
        -b "$ref_dir/TAIR10_1.CHG.bed" > \
        "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.CHG.bed"
    bedtools intersect -sorted -a "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.pileup.filt.bed" \
        -b "$ref_dir/TAIR10_1.CHH.bed" > \
        "$bam_dir/sorted_SQK-NBD114-96_barcode${barcode}.meth.tair10.CHH.bed"

    echo "Finished processing barcode $barcode."
    echo "----------------------------------------"
done

echo "All barcodes processed successfully!"