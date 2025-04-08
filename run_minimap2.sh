#!/bin/bash

# Directory containing the .sam files
directory="."

# Loop through each .fastq file in the directory
for fastq_file in "$directory"/col_2_pseudo/*.bam; do
	  # Extract the base name of the file (without extension)
	  #  base_name=$(basename "$fastq_file" .trim.fastq)
	  base_name=$(basename "$fastq_file" .bam)

	      # Define the output file names
	      sam_file="$directory/align_to_pseudo/${base_name}.sam"
	      bam_file="$directory/col_2_pseudo/${base_name}.bam"
		  sorted_bam_file="$directory/col_2_pseudo/sorted_${base_name}.bam"
          pseudoref_file="$directory/pseudoreference/TAIR10_1_$base_name.fasta"

		    # Run minimap2
		    echo -e "Running minimap2 on ${base_name} and ${pseudoref_file}"
		    minimap2 -ax map-ont -t 100 ${pseudoref_file} ${fastq_file} > ${sam_file}

		      # Run the samtools commands
		        samtools view -Sb "$sam_file" > "$bam_file"
			    samtools sort -o "$sorted_bam_file" -@ 100 "$bam_file"
			    samtools index -@ 100 "$sorted_bam_file"
            done