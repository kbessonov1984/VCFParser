#!/usr/bin/env python
import pandas as pd
import argparse, warnings
import os, re

#constants
support_ext = ["txt","tsv","vcf"]

def check_file_existance_and_type(path):

    abs_path = os.path.abspath(path)

    if not any([abs_path.endswith(ext) for ext in support_ext]):
        raise Exception("Unsupported file extension for {}. Supported file extensions are {}".format(path,support_ext))

    if os.path.exists(abs_path):
        return abs_path
    else:
        raise argparse.ArgumentTypeError("The {} can not be located".format(path))



def find_vcf_headerline_source(infile):
    source_vcf = ""
    for line in enumerate(open(file=infile).readlines()):
        regexmatch = re.search(r'source=(.+)', line[1])
        if regexmatch:
            source_vcf = regexmatch.group(1)
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

def get_input_type(path):
    file_extension = os.path.splitext(path)[1][1:].lower()
    return file_extension

def convert_tsv2vcf(tsvfilepath):
    vcf_list_for_df = [];
    vcf_df_header = ["CHROM","POS","ID","REF","ALT","TYPE","QUAL","FILTER", "INFO","FORMAT","SAMPLE"]
    with open(tsvfilepath) as fp:
        tsvheader = fp.readline().split()
        for line in fp:
            line = line.split()
            CHROM = line[0]
            POS = int(line[1])
            ID = '.'
            REF = line[2]
            ALT = line[3]

            var_type = 'SNP'
            if ALT[0] == '+':
                ALT = REF + ALT[1:]
                var_type = 'INS'
            elif ALT[0] == '-':
                REF += ALT[1:]
                ALT = line[2]
                var_type = 'DEL'
            QUAL = '.'
            pass_test = line[13]
            if pass_test == 'TRUE':
                FILTER = 'PASS'
            else:
                FILTER = 'FAIL'
            INFO = 'DP=' + line[11]
            FORMAT = 'GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ'
            SAMPLE = '1:' + line[4] + ':' + line[5] + ':' + line[6] + ':' + line[7] + ':' + line[8] + ':' + line[
                9] + ':' + line[10]
            vcf_list_for_df.append([CHROM,POS,ID,REF,ALT,var_type,QUAL,FILTER,INFO,FORMAT,SAMPLE])


    if 'REGION' and 'POS' not in tsvheader:
        raise Exception("Invalid TSV header input file {}".format(tsvfilepath))

    return pd.DataFrame(vcf_list_for_df,columns=vcf_df_header)

if __name__ == '__main__':

    parser = argparse.ArgumentParser("VCFparser.py parses iVar (https://andersen-lab.github.io/ivar/html/manualpage.html) "
                                     "TSV or VCF output files and filters out snvs linked to by the VOC\n")
    parser.add_argument('-i', '--input', required=True,
                        type=check_file_existance_and_type,
						help="List of ivar_variants.vcf or ivar_variants.tsv to summarise")
    parser.add_argument('-voc', '--voc_name',
						 required=True, help="Variant of Concern name (UK, SA) "
                        )
    parser.add_argument('-r', '--ref_meta', required=False, default='data/cov_lineage_variants.tsv',
                        type=check_file_existance_and_type,
						help="Path to metadata TSV file containing info on the key mutations")
    parser.add_argument('--signature_snvs_only', required=False, action='store_true',
                        help="Check VCF for only signature/official snvs linked to a VOC")
    parser.add_argument('--stat_filter_snvs', required=False, action='store_true',
                        help="Filter snvs based on statistical significance")

    args = parser.parse_args()

    inputtype = get_input_type(args.input)

    if inputtype == "vcf":
        skip_line_n, vcf_source = find_vcf_headerline_source(args.input)
        print("VCF source {}\nSelected VOC:{}".format(vcf_source, args.voc_name))
        vcf_df = pd.read_csv(args.input, sep="\t", skiprows=skip_line_n)
    elif inputtype == "tsv":
        vcf_df = convert_tsv2vcf(args.input)
        vcf_source = "TSV File"
    else:
        raise Exception("Unsupported input type")
    vcf_df.to_csv("test.tsv", sep="\t",index=False)

    if vcf_df.empty:
        raise  Exception("Empty input vcf_df DataFrame. Input parsing failed")


    input_file_name = os.path.basename(args.input)
    output_file_name = str(input_file_name.split(".vcf")[0])+"."+args.voc_name+".trimmed.vcf"

    VOCmeta_df = pd.read_csv(args.ref_meta, sep="\t")
    VOCmeta_df = VOCmeta_df[VOCmeta_df.VOC == args.voc_name]

    if VOCmeta_df.empty:
        raise Exception("Selected VOC {}  is not available in metadata file {}".format(args.voc_name,args.ref_meta))

    if vcf_source == 'iVar':
        VOCmeta_df.loc[VOCmeta_df["Type"] == "Del", "Position"] += 1


    if args.stat_filter_snvs:
        VOCmeta_df = VOCmeta_df[VOCmeta_df["SignatureSNV"] == True]

    if not type(vcf_df.loc[0,"POS"]) == type(VOCmeta_df.loc[0,"Position"]):
        raise Exception("vcf_df[\"POS\"] and VOCmeta_df[\"Position\"] types do not match. Check var type conversions")

    VOCmetaNotFound = VOCmeta_df[~VOCmeta_df["Position"].isin(vcf_df["POS"])]
    print("{} SNVs not found:\n {}".format(VOCmetaNotFound.shape[0],
        VOCmetaNotFound[["NucName","AAName","Position"]].to_string(index=False)))

    #filter 1: by posistion
    vcf_selected_idx = vcf_df["POS"].isin(VOCmeta_df["Position"]).to_list()
    vcf_df = vcf_df[vcf_selected_idx]
    if len(vcf_df.index) == 0:
        warnings.warn("Zero SNVs found in VCF. Might be an interesting sample or issue with input ... ")


    if args.stat_filter_snvs and "FILTER" in vcf_df.columns:
        vcf_df = vcf_df[vcf_df["FILTER"] == "PASS"]


    vcf_df_cleaned = remove_duplicated_vcf_snvs(vcf_df,VOCmeta_df)
    if vcf_df_cleaned.shape[0] == 0:
        vcf_df_cleaned.loc[0,"CHROM"] = "NO MATCHING SNVs"

    print("Writing out {} snvs to VCF".format(vcf_df_cleaned.shape[0]))
    vcf_df_cleaned.to_csv(output_file_name,sep="\t",index=False, mode="w")
    print("Trimmed VCF with VOC snvs is written to \"{}\"".format(output_file_name))
    print("Writting SNVs not found metadata tsv file \"{}\"".format("snvs_not_found.tsv"))
    VOCmetaNotFound.to_csv("snvs_not_found.tsv",sep="\t",index=False)
    print("Done")

