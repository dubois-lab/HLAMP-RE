#!/usr/bin/env python3

import csv
import gzip
import os
import sys
from collections import Counter

# ============================================================
# Configuration
# ============================================================

SAMPLE_SHEET = "./samplesheet_blank.csv"

VEP_DIR = "./output/vep"

OUTPUT_DIR = "./output/cbioportal"
GENOMIC_OUT = os.path.join(OUTPUT_DIR, "genomic_data.tsv")
CLINICAL_OUT = os.path.join(OUTPUT_DIR, "clinical_data.tsv")

CANCER_TYPE = "LUSC"

# Required cBioPortal mutation spectrum order
SPECTRUM_CLASSES = [
    "C>A",
    "C>G",
    "C>T",
    "T>A",
    "T>C",
    "T>G",
]

# VEP consequence -> cBioPortal alteration type
CONSEQUENCE_TO_TYPE = {
    "missense_variant": "MISSENSE",
    "stop_gained": "TRUNC",
    "frameshift_variant": "TRUNC",
    "stop_lost": "TRUNC",
    "start_lost": "TRUNC",
    "splice_donor_variant": "SPLICE",
    "splice_acceptor_variant": "SPLICE",
    "inframe_insertion": "INFRAME",
    "inframe_deletion": "INFRAME",
    "synonymous_variant": "SILENT",
}

# Lower number = higher priority
CONSEQUENCE_PRIORITY = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "missense_variant": 8,
    "inframe_insertion": 9,
    "inframe_deletion": 10,
    "synonymous_variant": 20,
}


# ============================================================
# Sample sheet
# ============================================================

def get_sample_ids(sample_sheet):
    """Read sample IDs from the sample sheet."""

    with open(sample_sheet, "r", newline="") as f:
        sample = f.read(4096)
        f.seek(0)

        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;")
        except csv.Error:
            dialect = csv.excel

        reader = csv.DictReader(f, dialect=dialect)

        if reader.fieldnames is None:
            raise ValueError(
                "Sample sheet does not contain a header."
            )

        possible_names = [
            "sampleID",
            "SampleID",
            "sample_id",
            "Sample",
            "sample",
        ]

        sample_column = None
        for name in possible_names:
            if name in reader.fieldnames:
                sample_column = name
                break

        if sample_column is None:
            raise ValueError(
                "Could not find a sample ID column.\n"
                f"Columns found: {reader.fieldnames}\n"
                "Expected one of: "
                + ", ".join(possible_names)
            )

        sample_ids = []
        for row in reader:
            sid = row[sample_column].strip()
            if sid:
                sample_ids.append(sid)

    return sample_ids


# ============================================================
# VEP parsing
# ============================================================

def parse_csq_format(header_line):
    """Extract VEP CSQ field names from the VCF header."""

    prefix = "Format: "
    if prefix not in header_line:
        raise ValueError(
            "Could not find VEP CSQ Format in VCF header."
        )

    fmt = header_line.split(prefix, 1)[1].strip()
    return fmt.split("|")


def parse_csq_entry(csq, csq_fields):
    """Convert one VEP CSQ annotation into a dictionary."""

    values = csq.split("|")

    if len(values) < len(csq_fields):
        values += [""] * (len(csq_fields) - len(values))

    return dict(zip(csq_fields, values))


def normalize_consequence(consequence):
    """
    VEP can report multiple consequences, e.g.
    missense_variant&splice_region_variant.
    Select the highest-priority recognized consequence.
    """

    consequences = consequence.split("&")

    recognized = [
        c for c in consequences
        if c in CONSEQUENCE_TO_TYPE
    ]

    if not recognized:
        return None

    return min(
        recognized,
        key=lambda c: CONSEQUENCE_PRIORITY.get(c, 999)
    )


def select_best_annotation(csq_string, csq_fields):
    """
    Select the most appropriate protein-coding VEP annotation.

    Preference:
      1. protein-coding transcript
      2. canonical transcript
      3. highest consequence priority
    """

    candidates = []

    for annotation in csq_string.split(","):
        if not annotation:
            continue

        csq = parse_csq_entry(annotation, csq_fields)

        if csq.get("BIOTYPE", "") != "protein_coding":
            continue

        consequence = normalize_consequence(
            csq.get("Consequence", "")
        )

        if consequence is None:
            continue

        canonical = (
            csq.get("CANONICAL", "").upper() == "YES"
        )

        priority = CONSEQUENCE_PRIORITY.get(
            consequence,
            999
        )

        candidates.append(
            (
                0 if canonical else 1,
                priority,
                csq
            )
        )

    if not candidates:
        return None

    candidates.sort(key=lambda x: (x[0], x[1]))
    return candidates[0][2]


# ============================================================
# Mutation spectrum
# ============================================================

def mutation_spectrum(ref, alt):
    """
    Convert an SNV into one of the six pyrimidine-based
    mutation classes used by cBioPortal.

    Returns None for non-SNVs.
    """

    ref = ref.upper()
    alt = alt.upper()

    if len(ref) != 1 or len(alt) != 1:
        return None

    if ref not in "ACGT" or alt not in "ACGT":
        return None

    if ref == alt:
        return None

    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }

    # Represent substitutions relative to the pyrimidine base.
    if ref in ("A", "G"):
        ref = complement[ref]
        alt = complement[alt]

    spectrum = f"{ref}>{alt}"

    if spectrum not in SPECTRUM_CLASSES:
        return None

    return spectrum


# ============================================================
# Protein alteration
# ============================================================

def hgvsp_to_alteration(hgvsp):
    """
    Example:
        ENSP00000369497.4:p.Val600Glu
    becomes:
        p.Val600Glu
    """

    if not hgvsp:
        return None

    if ":" in hgvsp:
        hgvsp = hgvsp.split(":", 1)[1]

    return hgvsp


# ============================================================
# Process one sample
# ============================================================

def process_sample(sid, genomic_writer):

    vcf_path = os.path.join(
        VEP_DIR,
        f"{sid}_tumor_only.vep.vcf.gz"
    )

    if not os.path.exists(vcf_path):
        print(
            f"WARNING: VCF not found: {vcf_path}",
            file=sys.stderr
        )
        return 0, Counter()

    mutation_count = 0
    spectrum_counts = Counter()
    csq_fields = None

    with gzip.open(vcf_path, "rt") as vcf:

        for line in vcf:

            if line.startswith("##INFO=<ID=CSQ"):
                if "Format: " in line:
                    csq_fields = parse_csq_format(line)
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            ref = fields[3]
            alt_field = fields[4]
            filt = fields[6]
            info = fields[7]

            if filt not in ("PASS", "."):
                continue

            csq_string = None
            for item in info.split(";"):
                if item.startswith("CSQ="):
                    csq_string = item[4:]
                    break

            if not csq_string:
                continue

            if csq_fields is None:
                raise ValueError(
                    f"No CSQ format found in {vcf_path}"
                )

            for alt in alt_field.split(","):

                annotation = select_best_annotation(
                    csq_string,
                    csq_fields
                )

                if annotation is None:
                    continue

                consequence = normalize_consequence(
                    annotation.get("Consequence", "")
                )

                if consequence is None:
                    continue

                alteration_type = CONSEQUENCE_TO_TYPE[
                    consequence
                ]

                # Exclude synonymous variants from the
                # main OncoPrint.
                if alteration_type == "SILENT":
                    continue

                gene = annotation.get(
                    "SYMBOL",
                    ""
                ).strip()

                # Fallback to Ensembl gene ID.
                if not gene:
                    gene = annotation.get(
                        "Gene",
                        ""
                    ).strip()

                if not gene:
                    continue

                alteration = hgvsp_to_alteration(
                    annotation.get("HGVSp", "")
                )

                if not alteration:
                    alteration = f"{ref}>{alt}"

                genomic_writer.writerow([
                    sid,
                    gene,
                    alteration,
                    alteration_type
                ])

                mutation_count += 1

                spectrum = mutation_spectrum(ref, alt)

                if spectrum:
                    spectrum_counts[spectrum] += 1

    return mutation_count, spectrum_counts


# ============================================================
# Main
# ============================================================

def main():

    print(f"Reading sample sheet: {SAMPLE_SHEET}")

    sample_ids = get_sample_ids(SAMPLE_SHEET)

    print(f"Found {len(sample_ids)} samples.")

    if not sample_ids:
        raise ValueError("No sample IDs were found.")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with open(
        GENOMIC_OUT,
        "w",
        newline=""
    ) as genomic_file:

        genomic_writer = csv.writer(
            genomic_file,
            delimiter="\t",
            lineterminator="\n"
        )

        genomic_writer.writerow([
            "sample",
            "HUGO_SYMBOL",
            "ALTERATION",
            "MUTATION_TYPE"
        ])

        clinical_rows = []

        for i, sid in enumerate(sample_ids, start=1):

            print(
                f"[{i}/{len(sample_ids)}] {sid}"
            )

            mutation_count, spectrum_counts = (
                process_sample(
                    sid,
                    genomic_writer
                )
            )

            # cBioPortal format:
            # C>A/C>G/C>T/T>A/T>C/T>G
            spectrum_string = "/".join(
                str(
                    spectrum_counts.get(
                        spectrum,
                        0
                    )
                )
                for spectrum in SPECTRUM_CLASSES
            )

            clinical_rows.append([
                sid,
                CANCER_TYPE,
                mutation_count,
                spectrum_string
            ])

            print(
                f"    mutations = {mutation_count}"
            )
            print(
                f"    spectrum  = {spectrum_string}"
            )

    with open(
        CLINICAL_OUT,
        "w",
        newline=""
    ) as clinical_file:

        clinical_writer = csv.writer(
            clinical_file,
            delimiter="\t",
            lineterminator="\n"
        )

        clinical_writer.writerow([
            "sample",
            "CANCER_TYPE",
            "MUTATION_COUNT(lognumber)",
            "MUTATION_SPECTRUM(C>A/C>G/C>T/T>A/T>C/T>G)"
        ])

        clinical_writer.writerows(clinical_rows)

    print()
    print("======================================")
    print("Finished")
    print("======================================")
    print(f"Genomic output : {GENOMIC_OUT}")
    print(f"Clinical output: {CLINICAL_OUT}")


if __name__ == "__main__":
    main()
