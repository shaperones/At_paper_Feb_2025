#!/bin/bash

# Directory containing the .sam files
directory="./"

# Loop through each .fastq file in the directory
for fastq_file in "$directory"samples/rep2/trim/*.trim.fastq; do
	  # Extract the base name of the file (without extension)
	  #  base_name=$(basename "$fastq_file" .trim.fastq)
	  base_name=$(basename "$fastq_file" .trim.fastq)

	      # Define the output file names
	      sam_file="$directory/${base_name}.sam"
	        bam_file="$directory/${base_name}.bam"
		  sorted_bam_file="$directory/sorted_${base_name}.bam"
		    nanotei_output="$directory/${base_name}.nanotei"

		    # Run minimap2
		    echo -e "Running minimap2 on ${base_name}"
		    minimap2 -ax map-ont -t 100 "./TAIR10_1.fna" ${fastq_file} > ${sam_file}

		      # Run the samtools commands
		        samtools view -Sb "$sam_file" > "$bam_file"
			  samtools sort -o "$sorted_bam_file" -@ 100 "$bam_file"
			    samtools index -@ 100 "$sorted_bam_file"

			      # Run the nanotei command
			        python3 ~/nanotei/nanotei.py "$sorted_bam_file" "$fastq_file" './TAIR10_1.fna' "./out_nanotei" "./Araprot11_TEs_fix.bed" "${base_name}.nanotei" --fpv --bed --minpvalue 0.05 -ovt 0.3
			done

