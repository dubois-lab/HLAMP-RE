import pysam
import csv
import argparse
import pandas as pd

def extract_het_info(input_vcf, output_file1):
    """
    Extract heterozygous variant information from a Clair3 VCF file.
    Parameters:
    -----------
    input_vcf : str
        Path to the input VCF file
    output_file1 : str
        Path to the output tab-delimited file
    """
    # Open the VCF file
    vcf_reader = pysam.VariantFile(input_vcf)
    # Open output file
    with open(output_file1, 'w') as out_f:
        # Write header
        out_f.write("seqnames\tstart\tend\talt.count\tref.count\n")
        # Iterate through VCF records
        for record in vcf_reader:
            # Check if the variant is heterozygous
            if record.samples[0]['GT'] != (1, 1):
                # Extract chromosome name
                seqnames = record.chrom
                # Start and end positions
                start = record.start
                end = start+1
                # Extract tumor and normal allele counts
                try:
                    alt_count = record.samples[0]['AD'][1]  # Alternate allele count for tumor
                    ref_count = record.samples[0]['AD'][0]  # Reference allele count for tumor
                    # Write to output file
                    out_f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(seqnames, start, end, alt_count, ref_count))
                
                except (IndexError, KeyError) as e:
                    print("Could not extract allele depths for record: {}. Error: {}".format(record, e))


def create_tumor_normal_hets(tumor_hets, normal_hets, output_file):
    # Load normal_hets file into a dictionary for fast lookup
    normal_data = {}
    with open(normal_hets, 'r') as f2:
        reader = csv.DictReader(f2, delimiter='\t')
        for row in reader:
            key = (row['seqnames'], row['start'], row['end'])
            normal_data[key] = (row['alt.count'], row['ref.count'])

    # Read tumor_hets file and write into a combined output_file
    with open(tumor_hets, 'r') as f1, open(output_file, 'w') as out:
        tumor_data = csv.DictReader(f1, delimiter='\t')
        fieldnames = ['seqnames', 'start', 'end', 'alt.count.t', 'ref.count.t', 'alt.count.n', 'ref.count.n']
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in tumor_data:
            key = (row['seqnames'], row['start'], row['end'])
            if key in normal_data:
                alt_count_n, ref_count_n = normal_data[key]
                writer.writerow({
                    'seqnames': row['seqnames'],
                    'start': row['start'],
                    'end': row['end'],
                    'alt.count.t': row['alt.count'],
                    'ref.count.t': row['ref.count'],
                    'alt.count.n': alt_count_n,
                    'ref.count.n': ref_count_n
                })
            else:
                alt_count_n = 0
                ref_count_n = 5 # median of normal read depth
                writer.writerow({
                    'seqnames': row['seqnames'],
                    'start': row['start'],
                    'end': row['end'],
                    'alt.count.t': row['alt.count'],
                    'ref.count.t': row['ref.count'],
                    'alt.count.n': alt_count_n,
                    'ref.count.n': ref_count_n
                })


def main():
    parser = argparse.ArgumentParser(description='Process VCF files for a specific tid.')
    parser.add_argument('tid', type=str, help='The tid to process')
    parser.add_argument('vcf', type=str, help='The tumor SNV vcf')
    args = parser.parse_args()

    print("TID: {}".format(args.tid))

    input_vcf_tumor = '{}'.format(args.vcf)
    input_pon_hets_counts = '../02_SNV_SV_calling/SNV_PoN_hets/PoN_hets_counts.tsv'

    output_file1_tumor = './inputdata/hets/intermediate_results/{}_tumor_hets_pysam.tsv'.format(args.tid)
    output_file3_combined = './inputdata/hets/intermediate_results/{}_hets_pysam.tsv'.format(args.tid)

    extract_het_info(input_vcf_tumor, output_file1_tumor)
    create_tumor_normal_hets(output_file1_tumor, input_pon_hets_counts, output_file3_combined)

if __name__ == "__main__":
    main()

