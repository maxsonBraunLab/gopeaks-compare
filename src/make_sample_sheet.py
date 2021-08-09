#!/usr/bin/env python3

import pandas as pd
import sys
import os
import re

def main():

    # get list of fastq files in raw dir
    raw = sys.argv[1]
    outfile = "samplesheet.tsv"

    samps = [os.path.join(raw,f) for f in os.listdir(raw)]
    regex = re.compile(".+\.fastq\.gz")
    sfilt = [s for s in samps if regex.search(s)]

    # make df 
    df=pd.DataFrame(sfilt, columns=["file"])
    df["sample"] = df["file"].apply(lambda x: "_".join(os.path.basename(x).split("_")[0:3]))
    df["mark"] = df["file"].apply(lambda x: x.split("_")[2])
    df["condition"] = df['sample'].apply(lambda x: x.split("_")[0])

    df["read"] = df["file"].apply(lambda x: "R1" if "R1" in x else "R2")
    df["igg"] = df["sample"]
    df["peak"] = df["mark"].apply(lambda x: "NA" if "IgG" in x else "narrow")

    df=df.pivot_table(index=['sample','igg','mark','condition', "peak"], columns='read', values='file',aggfunc=sum).reset_index()
    df=df[["sample","R1","R2","mark","condition","igg", "peak"]]
    df.columns.name = None

    # write sample sheet to file
    df.to_csv(outfile, sep='\t', index=False)

    # write deseq2 metadata to file
    dmet=df.copy()
    dmet=dmet[['sample','condition']]
    dmet.loc[:,('sample')] = dmet.loc[:,('sample')].apply(lambda x: x.split("_")[0])
    dmet.drop_duplicates(inplace=True)
    dmet.to_csv('src/deseq2_metadata.csv', index=False)

if __name__ == "__main__": 
    main()


# changes
# assume input is {condition}_{replicate}_{mark}_{R1|R2}.fastq.gz
# add a 'peak' column for macs2 broad or narrow option.