# VCFParser

Extract non-duplicated SNVs from a VCF file associated with a given SARS-COV2 Variable of Concern (VOC)
Append new SNVs to the `data/cov_lineage_variants.tsv` file

### Requirements
* pandas
* python >= 3
* matplotlib >= 3 

### Usage

```
$ python vcfparser.py -h
usage: VCFparser.py parses iVar (https://andersen-lab.github.io/ivar/html/manualpage.html) TSV or VCF output files and filters out snvs linked to by the VOC

       [-h] -i INPUT [INPUT ...] [-bam BAM_FILES [BAM_FILES ...]] -voc
       VOC_NAMES [-r REF_META] [--signature_snvs_only] [--key_snvs_only]
       [--stat_filter_snvs] --subplots_mode SUBPLOTS_MODE
       [--min_snv_freq_threshold [0-1]] [--annotate]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        List of ivar_variants.vcf or ivar_variants.tsv to
                        summarise
  -bam BAM_FILES [BAM_FILES ...], --bam_files BAM_FILES [BAM_FILES ...]
                        Optionally provide a list of corresponding bam files
                        in THE SAME ORDER asfiles provided for the -i
                        parameter
  -voc VOC_NAMES, --voc_names VOC_NAMES
                        List of Variants of Concern names (e.g. UK, SA,
                        Brazil, Nigeria)
  -r REF_META, --ref_meta REF_META
                        Path to metadata TSV file containing info on the key
                        mutations
  --signature_snvs_only
                        Check VCF for only signature/official snvs linked to a
                        VOC
  --key_snvs_only       Check VCF for only the key (S-gene associated) snvs
                        linked to a VOC
  --stat_filter_snvs    Filter snvs based on statistical significance (i.e. QC
                        PASS flag)
  --subplots_mode SUBPLOTS_MODE
                        How to plot multiple plots (onerow, onecolumn,
                        oneplotperfile)
  --min_snv_freq_threshold [0-1]
                        Set minimum SNV frequency threshold to display
                        (default: 0)
  --annotate            Annotate heatmap with SNV frequency values
```
