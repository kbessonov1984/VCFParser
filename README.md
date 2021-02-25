# VCFParser

Extract non-duplicated SNVs from a VCF file associated with a given SARS-COV2 Variable of Concern (VOC)
Append new SNVs to the `data/cov_lineage_variants.tsv` file

### Requirements
* pandas
* python >= 3


### Usage

```
$ python vcfparser.py -h
usage: VCFparser.py parses iVar (https://andersen-lab.github.io/ivar/html/manualpage.html) TSV or VCF output files and filters out snvs linked to by the VOC

       [-h] -i INPUT [INPUT ...] -voc VOC_NAMES [-r REF_META]
       [--signature_snvs_only] [--stat_filter_snvs]
       [--subplots_orientation SUBPLOTS_ORIENTATION]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        List of ivar_variants.vcf or ivar_variants.tsv to
                        summarise
  -voc VOC_NAMES, --voc_names VOC_NAMES
                        List of Variants of Concern names (e.g. UK, SA,
                        Brazil, Nigeria)
  -r REF_META, --ref_meta REF_META
                        Path to metadata TSV file containing info on the key
                        mutations
  --signature_snvs_only
                        Check VCF for only signature/official snvs linked to a
                        VOC
  --stat_filter_snvs    Filter snvs based on statistical significance
  --subplots_orientation SUBPLOTS_ORIENTATION
                        List of Variants of Concern names (e.g. UK, SA,
                        Brazil, Nigeria)
```
