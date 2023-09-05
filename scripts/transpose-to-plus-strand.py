import os
import sys
import shutil 
import subprocess
import tempfile
import logging
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
import argparse

# initialize logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# read in original transcripts 
def read_fasta_file(fasta_file):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    logger.info(f"read {len(fasta_dict)} records from {fasta_file}")
    return fasta_dict

# logic for creating gtf file and transcript file -- used by create_transcripts 
def write_transcripts_and_exons(chrom, seq_record, strand, gtf_out, fasta_out, fasta_strand_out):
    chrom_length = len(seq_record)
    transcript_line = f"{chrom}\t.\ttranscript\t1\t{chrom_length}\t.\t{strand}\t.\tgene_id \"{chrom}_{strand}\"; gene_name \"{chrom}_{strand}\"; transcript_id \"{chrom}_{strand}\"; transcript_biotype \"chromosome_strand\";\n"
    exon_line = f"{chrom}\t.\texon\t1\t{chrom_length}\t.\t{strand}\t.\tgene_id \"{chrom}_{strand}\"; transcript_id \"{chrom}_{strand}\";\n"
    gtf_out.write(transcript_line + exon_line)
    seq = seq_record.seq if strand == '+' else seq_record.seq.reverse_complement()
    fasta_out.write(f">{chrom}_{strand}\n{str(seq)}\n")
    fasta_strand_out.write(f">{chrom}_{strand}\n{str(seq)}\n")

# create plus-strand reference fasta and gtf file 
def create_transcripts(fasta_dict, output_gtf, output_fasta, output_fasta_plus, output_fasta_minus):
    with open(output_gtf, "w") as gtf_out, open(output_fasta, "w") as fasta_out, \
         open(output_fasta_plus, "w") as fasta_plus_out, open(output_fasta_minus, "w") as fasta_minus_out:
        for chrom, seq_record in fasta_dict.items():
            write_transcripts_and_exons(chrom, seq_record, '+', gtf_out, fasta_out, fasta_plus_out)
            write_transcripts_and_exons(chrom, seq_record, '-', gtf_out, fasta_out, fasta_minus_out)
    logger.info(f"written transcripts and exons to {output_gtf}, {output_fasta}, {output_fasta_plus}, and {output_fasta_minus}")

# filter for priamry, mapped reads, and split these sequences by the strand to which they were aligned
def split_bam_by_strand(input_bam, plus_bam, minus_bam):

    plus_count = 0
    minus_count = 0

    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
         pysam.AlignmentFile(plus_bam, "wb", header=in_bam.header) as plus_out, \
         pysam.AlignmentFile(minus_bam, "wb", header=in_bam.header) as minus_out:
        for read in in_bam:
            if not read.is_secondary and not read.is_supplementary:
                if read.is_reverse:
                    minus_out.write(read)
                    minus_count += 1
                else:
                    plus_out.write(read)
                    plus_count += 1

    logger.info(f"Split {plus_count} reads to {plus_bam} and {minus_count} reads to {minus_bam}")

    # sort and index the BAM files
    for bam_file in [plus_bam, minus_bam]:
        sorted_bam = bam_file.replace('.bam', '.sorted.bam')
        pysam.sort('-o', sorted_bam, '-@', '104', bam_file)
        
        # move the sorted file back to the original bam_file
        os.rename(bam_file, bam_file.replace('.bam', '.original.bam'))  # rename original file
        shutil.move(sorted_bam, bam_file)  # move sorted file to original bam_file name
        pysam.index(bam_file)
        logger.info(f"Sorted and indexed {bam_file}")

# extract the original read sequences from the bam file 
def extract_sequences_from_bam(input_bam, output_fastq):
    
    count = 0
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, open(output_fastq, "w") as fastq_out:
        for read in in_bam:
            seq = read.query_sequence
            qual = ''.join([chr(q + 33) for q in read.query_qualities])

            # reverse complement the sequence and reverse quality if the read is from the minus strand
            # because
            # the sequence in the BAM file is reverse complimented (compared to the original fastq) if the alignment is on the minus strand 
            # see: https://www.biostars.org/p/198051/
            if read.is_reverse:
                seq = Seq(seq).reverse_complement()  # rev-comp minus strand 
                qual = qual[::-1]

            fastq_out.write(f"@{read.query_name}\n{str(seq)}\n+\n{qual}\n")
            count += 1
    logger.info(f"Written {count} sequences to {output_fastq}")


# align the plus-stranded reads to the plus-strand reference, and the minus-strand reads to the transposed minus strand reference 
def minimap_align(input_sequence, plus_strand_genome, output_bam, minimap_params):

    # prepare minimap2 args, either default or user-supplied 
    minimap_args = minimap_params.split()

    # Prepare the command for minimap2 alignment
    minimap_cmd = ["minimap2"] + minimap_args + [plus_strand_genome, input_sequence]

    # redirect the alignment to the output_bam
    with open(output_bam, "w") as out_file:
        try:
            subprocess.run(minimap_cmd, stdout=out_file, check=True)
        except subprocess.CalledProcessError as e:
            logger.error("minimap2 failed with message: " + str(e))
            sys.exit(1)

    logger.info(f"Aligned sequences in {input_sequence} to {plus_strand_genome}, output in {output_bam}")

# merge alignments to the plus strand and transposed minus strand 
def merge_bam_files(input_bam1, input_bam2, output_bam):
    merge_cmd = ["samtools", "merge", "-f", output_bam, input_bam1, input_bam2]
    try:
        subprocess.run(merge_cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("samtools merge failed with message: " + str(e))
        sys.exit(1)
    logger.info(f"Merged {input_bam1} and {input_bam2} into {output_bam}")

# process by strand
def process_strand(input_bam, plus_strand_ref_fasta, minus_strand_ref_fasta, output_dir, minimap_params):

    # temporary BAM file for plus strand
    plus_strand_bam = tempfile.NamedTemporaryFile(dir=output_dir, delete=False, suffix=".bam").name

    # temporary BAM file for minus strand
    minus_strand_bam = tempfile.NamedTemporaryFile(dir=output_dir, delete=False, suffix=".bam").name

    # split reads by strand
    split_bam_by_strand(input_bam, plus_strand_bam, minus_strand_bam)

    aligned_bams = []  # List to store the aligned BAM files

    # process plus and minus strand separately 
    for strand_bam, ref_fasta, strand in [(plus_strand_bam, plus_strand_ref_fasta, "+"), (minus_strand_bam, minus_strand_ref_fasta, "-")]:
        
        # create a temporary fastq file and write reads from the stranded file to the fastq 
        output_fastq = tempfile.NamedTemporaryFile(dir=output_dir, delete=False, suffix=".fastq").name
        extract_sequences_from_bam(strand_bam, output_fastq)

        # create a temporary BAM file for the aligned sequences and align sequences with minimap2
        aligned_bam = tempfile.NamedTemporaryFile(dir=output_dir, delete=False, suffix=".bam").name
        minimap_align(output_fastq, ref_fasta, aligned_bam, minimap_params)

        # append the new alignments to "aligned bams"
        aligned_bams.append(aligned_bam)  

    return aligned_bams  

def main():

    # get the directory path from output_merged_bam
    output_dir = os.path.dirname(args.output_merged_bam)

    # fetch the basename for the output fasta 
    base_fasta_name = os.path.splitext(args.output_fasta)[0]  

    try:
        # read FASTA file and create a dict of sequences 
        fasta_dict = read_fasta_file(args.input_fasta)

        # create a plus strand reference, including an output GTF and output FASTA 
        plus_strand_ref_fasta = base_fasta_name + "_plus.fasta"
        minus_strand_ref_fasta = base_fasta_name + "_minus.fasta"
        create_transcripts(fasta_dict, args.output_gtf, args.output_fasta, plus_strand_ref_fasta, minus_strand_ref_fasta)

        # run the process strand function 
        plus_aligned_bam, minus_aligned_bam = process_strand(args.input_bam, plus_strand_ref_fasta, minus_strand_ref_fasta, output_dir, args.minimap_params)

        # merge the two bam files of primary alignments
        merge_bam_files(plus_aligned_bam, minus_aligned_bam, args.output_merged_bam)

        # clean up temporary files 
        os.remove(plus_aligned_bam)
        os.remove(minus_aligned_bam)

    except Exception as e:
        logger.error(f"An error occurred: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('input_bam', type=str, help='Input BAM file')
    parser.add_argument('input_fasta', type=str, help='Input FASTA file')
    parser.add_argument('output_fasta', type=str, help='Output FASTA file')
    parser.add_argument('output_gtf', type=str, help='Output GTF file')
    parser.add_argument('output_merged_bam', type=str, help='Output merged BAM file')
    parser.add_argument('--minimap_params', type=str, default='-a -x splice -k 15 --for-only -t 104', help='Parameters for minimap alignment')

    args = parser.parse_args()
    main()
