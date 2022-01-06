# common holds pyhton function to be used in the snakefile
import pandas as pd
import yaml

# map samples to fastqs
def get_samples():
    """
    return list of samples from samplesheet.tsv
    """
    return list(st.index)

def get_marks():
    """
    return list of marks from samplesheet.tsv
    """
    return list(st['mark'])

def get_bowtie2_input(wildcards):
    """
    returns reads associated with a sample
    """
    r1=st.loc[wildcards.sample]['R1']
    r2=st.loc[wildcards.sample]['R2']
    return r1,r2

def get_reads():
    """
    get list of all reads
    """
    rlist=list(st['R1'])+list(st['R2'])
    rlist=[os.path.basename(f).split('.')[0] for f in rlist]
    return rlist

def get_igg(wildcards):
    """
    Returns the igg file for the sample unless
    the sample is IgG then no control file is used.
    """
    if config['USEIGG']:
        igg=st.loc[wildcards.sample]['igg']
        iggbam=f'data/ban/{igg}.ban.sorted.markd.bam'
        isigg="IgG" in wildcards.sample
        if not isigg:
            return f'{iggbam}'
        else:
            # return ""
            return f'{iggbam}'
    else:
        return ""

def macs2_igg(wildcards):
    """
    Returns the igg file for the sample unless
    the sample is IgG then no control file is used.
    """ 
    if config['USEIGG']:
        igg=st.loc[wildcards.sample]['igg']
        iggbam=f'data/ban/{igg}.ban.sorted.markd.bam'
        isigg="IgG" in wildcards.sample
        if not isigg:
            return f'-c {iggbam}'
        else:
            return ""
    else:
        return ""

def macs2_peak(wildcards):
    """
    Returns the peak type (broad or narrow) for a sample.
    """
    peak=st.loc[wildcards.sample]['peak']
    if peak == 'broad':
        return "--broad"
    if peak == 'narrow':
        return ""

def seacr_igg(wildcards):
    """
    Returns the igg file for the sample unless
    the sample is IgG then no control file is used.
    """
    if config['USEIGG']:
        igg=st.loc[wildcards.sample]['igg']
        iggbam=f'data/seacr/bedgraph/{igg}.bedgraph'
        isigg="IgG" in wildcards.sample
        if not isigg:
            return f'{iggbam}'
        else:
            return config['SEACR_THRESHOLD']
    else:
        return ""

def seacr_norm(wildcards):
    """
    Returns "norm" if using SEACR with treatment + control situation
    Returns "non" if using SEACR with no treatment. E.g. just IgG
    """
    if config["IGG"] in wildcards.sample:
        return "non"
    else:
        return "norm"

def get_callpeaks(wildcards):
    """
    Returns the callpeaks input files
    """
    bam=f"data/ban/{wildcards.sample}.ban.sorted.markd.bam"
    bai=f"data/ban/{wildcards.sample}.ban.sorted.markd.bam.bai"
    # cp="scripts/gopeaks"
    return [bam,bai]

def gopeaks_igg(wildcards):
    """
    Returns the igg file for the sample unless
    the sample is IgG then no control file is used.
    """
    if config['USEIGG']:
        igg=st.loc[wildcards.sample]['igg']
        iggbam=f'data/ban/{igg}.ban.sorted.markd.bam'
        isigg="IgG" in wildcards.sample
        if not isigg:
            return f'-control {iggbam}'
        else:
            return ""
    else:
        return ""

def dynamic_range_input():
    """
    Input: list of all consensus peak files
    Method: Loop through each consensus file. Define method,condition,mark.
    Method: Use method,condition,mark to glob relevant bam file pairs.
    Output: One consensus file with one BAM file. Their condition,mark must match,
    Output: but must output all methods for a sample.
    """
    all_consensus = glob.glob("data/consensus/*.bed")
    consensus_file = os.path.basename(consensus).replace(".bed", "")
    method=consensus_file.split("_")[0]
    condition=consensus_file.split("_")[1]
    mark=consensus_file.split("_")[2]
    input_bams = glob.glob("data/ban/{condition}*{mark}*.bam".format(condition = condition, mark = mark))
    for bam in input_bams:
        return(consensus_file, bam)

def get_standard(wildcards):
    """
    Input: wildcard of a sample
    Method: use wildcard of a sample to grab standard in config.yml
    Output: the gold standard file
    """
    with open("src/config.yml", "r") as fi:
        cfg = yaml.safe_load(fi)
    return( cfg["STANDARDS"][wildcards.sample] )

def get_minwidth(wildcards):
    """
    Input: wildcard of a sample
    Output: change minwidth for H3K27Ac in Kasumi cells
    Note: only appropriate for this analysis
    """
    condition = str(wildcards.sample).split("_")[0]
    replicate = str(wildcards.sample).split("_")[1]
    mark = str(wildcards.sample).split("_")[2]
    if ((condition == "Kasumi") and (mark == "H3K27Ac")):
        return ""
    else:
        return ""

def get_groups():
    """
    Get {condition}_{mark} information without replicates using the samplesheet.
    """
    outdf = pd.DataFrame()
    cond_mark = st[['condition', 'mark']][st.mark != "IgG"].drop_duplicates().reset_index(drop = True) # condition_mark w/out replicates
    for method in all_methods:
        tmpdf = cond_mark.copy()
        tmpdf['method'] = method
        outdf = pd.concat([outdf, tmpdf])
    outdf = outdf.reset_index(drop = True)
    return(outdf)

def group_reps(wildcards):
    """
    Input: {condition}_{mark} groups
    Output: peak files with {method}_{condition}_{replicate}_{mark}
    """
    samples = list( st[(st.condition == wildcards.condition) & (st.mark == wildcards.mark)]['sample'] )
    if wildcards.method == "gopeaks":
        return([ "data/gopeaks/{s}.bed".format(s = i) for i in samples ])
    if wildcards.method == "macs2":
        return([ "data/macs2/{s}_peaks.narrowPeak".format(s = i) for i in samples ])
    if wildcards.method == "seacr-relaxed":
        return([ "data/seacr/{s}.relaxed.bed".format(s = i) for i in samples ])
    if wildcards.method == "seacr-stringent":
        return([ "data/seacr/{s}.stringent.bed".format(s = i) for i in samples ])

