#!/usr/bin/env python

__author__ = "Tiurin K."
__author_email__ = "tiurin.kn@gmail.com"

import os
import argparse
import sys
import pysam
import re
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import warnings

def get_args():
    ###get arguments from command line
    desc = (
        """"""
    )
    epi = """"""
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument("BAM", action="store", help='path to sorted .bam file from WGS')
    parser.add_argument("bed", action="store", help='path to .bed file with insertions positions')
    parser.add_argument("output_dir", action="store", help='path to output directory')
    parser.add_argument("prefix", action="store", help='output file prefix')
    parser.add_argument("--cpu", action="store", help='CPU number [default 1]')
    

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()

def extract_soft_clipped(seq, cigar, ref_start):
    """
    Extracts soft-clipped sequences from the read.
    """
    soft_clipped_seqs = []
    ref_pos = ref_start
    query_pos = 0
    
    # Parse CIGAR string
    for op, length in cigar:
        if op == 0:  # Match (M)
            ref_pos += length
            query_pos += length
        elif op == 4:  # Soft-clipping (S)
            soft_clipped_seq = seq[query_pos:query_pos + length]
            soft_clipped_seqs.append((ref_pos, soft_clipped_seq))
            query_pos += length
        else:
            query_pos += length if op in {1, 4} else 0
            ref_pos += length if op in {0, 2, 3} else 0
    
    return soft_clipped_seqs

def extract_clipped_part_from_reads(BAM, chr, start, end, out_fasta, ins_reads_lengths_aver):
    """
    Extracts and writes soft-clipped sequences within the given range to a FASTA file.
    """
    with open(out_fasta, 'w') as new, pysam.AlignmentFile(BAM, "rb") as samfile:
        for read in samfile.fetch(chr, start, end):
            if not read.cigartuples:  # Skip reads without a CIGAR string
                continue
            if not read.query_sequence:  # Skip reads without a CIGAR string
                continue
            seq = read.query_sequence
            cigar = read.cigartuples  # Parsed CIGAR
            read_name = read.query_name
            
            soft_clipped_parts = extract_soft_clipped(seq, cigar, read.reference_start)
            
            for clip_pos, soft_clipped_seq in soft_clipped_parts:
                if ins_reads_lengths_aver > 0:
                    if start <= clip_pos < end and len(soft_clipped_seq) > 100 and len(soft_clipped_seq) < ins_reads_lengths_aver + 100:
                        new.write(f'>{read_name}\n{soft_clipped_seq}\n')
                else:
                    if start <= clip_pos < end and len(soft_clipped_seq) > 100:
                        new.write(f'>{read_name}\n{soft_clipped_seq}\n')

def extract_insertions(seq, cigar, ref_start):
    insertion_seqs = []
    ref_pos = ref_start
    query_pos = 0
    
    # Parse CIGAR string
    for op, length in cigar:
        if op == 0:  # Match (M) or equal (E) => reference position moves forward
            ref_pos += length
            query_pos += length
        elif op == 1:  # Insertion (I) => insertion in the query sequence
            insertion_seq = seq[query_pos:query_pos + length]  # Get insertion sequence from query
            insertion_pos = ref_pos  # Insertion happens at the current reference position
            insertion_seqs.append((insertion_pos, insertion_seq))
            query_pos += length
        elif op == 2:  # Deletion (D) => no change in the query
            ref_pos += length
    
    return insertion_seqs

def extract_insertions_from_reads(BAM, chr, start, end, out_fasta):
    # Open the BAM file
    samfile = pysam.AlignmentFile(BAM, "rb")

    ins_reads_lengths = []
    
    with open(out_fasta, 'w') as new:
        for read in samfile.fetch(chr, start, end):
            if not read.cigartuples:  # Skip reads without a CIGAR string
                continue
            if not read.query_sequence:  # Skip reads without a CIGAR string
                continue
            seq = read.query_sequence
            cigar = read.cigar
            read_name = read.query_name
            
            insertions = extract_insertions(seq, cigar, read.reference_start)
                
            for insertion_pos, insertion_seq in insertions:
                if start <= insertion_pos < end and len(insertion_seq) > 100:
                    ins_reads_lengths.append(len(insertion_seq))
                    new_line = f'>{read_name}\n{insertion_seq}\n'
                    new.write(new_line)

    if len(ins_reads_lengths) != 0:
        ins_reads_lengths_aver = sum(ins_reads_lengths) / len(ins_reads_lengths)
    else:
        ins_reads_lengths_aver = 0
    
    samfile.close()

    return ins_reads_lengths_aver

def run_mafft(cpu, fasta_I, fasta_S, outpath, name_for_insertion, out_insertion):
    os.system('mkdir {0}'.format(outpath))
    if fasta_I == 'NO':
        os.system('mafft --quiet --thread {0} {1} > {2}/{3}.aln.fasta'.format(cpu, fasta_S, outpath, name_for_insertion))
    if fasta_S == 'NO':
        os.system('mafft --quiet --thread {0} {1} > {2}/{3}.aln.fasta'.format(cpu, fasta_I, outpath, name_for_insertion))
    if fasta_I != 'NO' and fasta_S != 'NO':
        os.system('cat {0} {1} > {2}/{3}_merg.fasta'.format(fasta_I, fasta_S, outpath, name_for_insertion))
        os.system('mafft --quiet --thread {0} {1}/{2}_merg.fasta > {1}/{2}.aln.fasta'.format(cpu, outpath, name_for_insertion))

    alignment = AlignIO.read(f'{outpath}/{name_for_insertion}.aln.fasta', "fasta")
    summary = AlignInfo.SummaryInfo(alignment)
    consensus = summary.dumb_consensus(0.5, "N")

    with open(out_insertion, 'w') as new:
        new.write(f">{name_for_insertion}\n{consensus}\n")

    os.system('rm {0}'.format(fasta_I))
    os.system('rm {0}'.format(fasta_S))
    os.system('rm -rf {0}'.format(outpath))

def main():
    print('Thanks for using xx')
    args = get_args()

    warnings.filterwarnings("ignore")

    if args.cpu:
        cpu = args.cpu
    else:
        cpu = 1

    with open(args.bed, 'r') as bed_file:
        for line in bed_file:
            line = line.strip().split('\t')
            chr, start, end, insertion_name = line[0], int(line[1]), int(line[2]), line[3]

            ins_reads_lengths_aver = extract_insertions_from_reads(args.BAM, chr, start - 5, end + 5, f'{args.output_dir}/{args.prefix}_{insertion_name}_I.fasta')
            extract_clipped_part_from_reads(args.BAM, chr, start - 2, end + 2, f'{args.output_dir}/{args.prefix}_{insertion_name}_S.fasta', ins_reads_lengths_aver)
            

            is_ins_reads = os.stat(f'{args.output_dir}/{args.prefix}_{insertion_name}_I.fasta').st_size == 0
            is_soft_clip_reads = os.stat(f'{args.output_dir}/{args.prefix}_{insertion_name}_S.fasta').st_size == 0

            if is_ins_reads == False and is_soft_clip_reads == False:
                run_mafft(cpu,
                         f'{args.output_dir}/{args.prefix}_{insertion_name}_I.fasta',
                         f'{args.output_dir}/{args.prefix}_{insertion_name}_S.fasta',
                         f'{args.output_dir}/tmp', insertion_name, f'{args.output_dir}/{args.prefix}_{insertion_name}_assembled.fasta')
            
            if is_ins_reads == False and is_soft_clip_reads == True:
                run_mafft(cpu,
                         f'{args.output_dir}/{args.prefix}_{insertion_name}_I.fasta',
                         'NO',
                         f'{args.output_dir}/tmp', insertion_name, f'{args.output_dir}/{args.prefix}_{insertion_name}_assembled.fasta')
                
            if is_ins_reads == True and is_soft_clip_reads == False:
                run_mafft(cpu,
                         'NO',
                         f'{args.output_dir}/{args.prefix}_{insertion_name}_S.fasta',
                         f'{args.output_dir}/tmp', insertion_name, f'{args.output_dir}/{args.prefix}_{insertion_name}_assembled.fasta')

    os.system('cat {0}/*_assembled.fasta > {0}/{1}_insertions.fasta'.format(args.output_dir, args.prefix))
    os.system('rm {0}/*_assembled.fasta'.format(args.output_dir))


if __name__ == "__main__":
    main()