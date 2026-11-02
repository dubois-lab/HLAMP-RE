import glob
import pandas as pd
import os

output_folder = "inputdata/coverage"
cnvkit_out = "../03_CNV_calling/inputs"
cnn_files = glob.glob(os.path.join(cnvkit_out, "*_tumor_5kbBin_UCSC.cnn"))

for sample in cnn_files:
    df = pd.read_csv(sample, sep="\t")
    # Calculate linear ratio from log2 ratio from cnvkit
    df['ratio'] = 2 ** df['log2']

    output_df = df[['chromosome', 'start', 'end', 'ratio']].copy()
    output_df.columns = ['seqnames', 'start', 'end', 'ratio']
    base_name = os.path.splitext(os.path.basename(sample))[0]
    output_csv = os.path.join(output_folder, base_name + "_cov.csv")

    output_df.to_csv(output_csv, index=False)
    print(f"Processed {sample} -> {output_csv}")

