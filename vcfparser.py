#!/usr/bin/env python
import pandas as pd
import argparse
import os, re


def check_file_existance(path):

    abs_path = os.path.abspath(path)
    if os.path.exists(abs_path):
        return abs_path
    else:
        raise argparse.ArgumentTypeError("The {} can not be located".format(path))

def find_vcf_headerline_source(infile):
    source_vcf = ""
    for line in enumerate(open(file=infile).readlines()):
        if re.findall(r'##source=iVar', line[1]):
            source_vcf = "iVar"
        if re.findall(r'CHROM',line[1]):
            return line[0], source_vcf


def remove_duplicated_vcf_snvs(vcf, meta):
    idx_dupl = vcf[vcf["POS"].duplicated(keep=False)].index.values.tolist()
    positions_dupl = set(vcf.loc[idx_dupl, "POS"].to_list())

    for position in positions_dupl:
        vcf_tmp = vcf[vcf["POS"] == position]
        Ref_allele = meta[meta["Position"]==position]["Ref"].values[0]
        Alt_allele = meta[meta["Position"] == position]["Alt"].values[0]
        idxs = vcf_tmp[(vcf_tmp["REF"] == Ref_allele) & (vcf_tmp["ALT"] == Alt_allele)].index.values


        for remove_idx in idxs:
            idx_dupl.remove(remove_idx)


    return (vcf.drop(idx_dupl))




if __name__ == '__main__':

    parser = argparse.ArgumentParser("VCFparser.py parses VCF file and outputs the defined "
                                     "key mutations in TSV or VCF\n")
    parser.add_argument('-i', '--input', required=True,
                        type=check_file_existance,
						help="List of ivar_variants.vcf to summarise")
    parser.add_argument('-voc', '--voc_name',
						 required=True, help="Variant of Concern name (UK, SA) "
                        )
    parser.add_argument('-r', '--ref_meta', required=False, default='data/cov_lineage_variants.tsv',
                        type=check_file_existance,
						help="Path to metadata TSV file containing info on the key mutations")
    parser.add_argument('--signature_snvs_only', required=False, action='store_true',
                        help="Check VCF for only signature/official snvs linked to VOCs")

    args = parser.parse_args()
    skip_line_n, vcf_source = find_vcf_headerline_source(args.input)
    print("VCF source {}\nSelected VOC:{}".format(vcf_source, args.voc_name))


    vcf_df = pd.read_csv(args.input, sep="\t", skiprows=skip_line_n)

    input_file_name = os.path.basename(args.input)
    output_file_name = str(input_file_name.split(".vcf")[0])+"."+args.voc_name+".trimmed.vcf"

    VOCmeta_df = pd.read_csv(args.ref_meta, sep="\t")
    VOCmeta_df = VOCmeta_df[VOCmeta_df.VOC == args.voc_name]

    if VOCmeta_df.empty:
        raise Exception("Selected VOC {}  is not available in metadata file {}".format(args.voc_name,args.ref_meta))

    if vcf_source == 'iVar':
        VOCmeta_df.loc[VOCmeta_df["Type"] == "Del", "Position"] += 1


    if args.signature_snvs_only:
        VOCmeta_df = VOCmeta_df[VOCmeta_df["SignatureSNV"] == True]


    #filter 1: by posistion
    vcf_selected_idx = vcf_df["POS"].isin(VOCmeta_df["Position"]).to_list()



    vcf_df = vcf_df[vcf_selected_idx]

    if "FILTER" in vcf_df.columns:
        vcf_df = vcf_df[vcf_df["FILTER"] == "PASS"]



    vcf_df_cleaned = remove_duplicated_vcf_snvs(vcf_df,VOCmeta_df)
    vcf_df_cleaned.to_csv(output_file_name,sep="\t",index=False, mode="w")
    print("Trimmed VCF with VOC snvs is written to {}".format(output_file_name))
    print("Done")

