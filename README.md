#VCFParser

# iVar tool variant report generator

Extract SNVs from a VCF file associated with a given SARS-COV2 Variable of Concern (VOC)
Append new SNVs to the `data/cov_lineage_variants.tsv` file

### Requirements
* pandas
* python >= 3


### Usage

```
$ python vcfparser.py -h
usage: VCFparser.py parses VCF file and outputs the defined key mutations in TSV or VCF

       [-h] -i INPUT -voc VOC_NAME [-r REF_META] [--signature_snvs_only]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        List of ivar_variants.vcf to summarise
  -voc VOC_NAME, --voc_name VOC_NAME
                        Variant of Concern name (UK, SA)
  -r REF_META, --ref_meta REF_META
                        Path to metadata TSV file containing info on the key
                        mutations
  --signature_snvs_only
                        Check VCF for only signature/official snvs linked to
                        VOCs
```
