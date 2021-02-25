#!/usr/bin/env python
import pandas as pd
import argparse, warnings
import os, re
import matplotlib.pyplot as plt

import  VOCheatmapper

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
    vcf_list_for_df = []
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

def csv_list(string):
   return string.split(',')

if __name__ == '__main__':
    parser = argparse.ArgumentParser("VCFparser.py parses iVar (https://andersen-lab.github.io/ivar/html/manualpage.html) "
                                     "TSV or VCF output files and filters out snvs linked to by the VOC\n")
    parser.add_argument('-i', '--input', required=True,
                        type=check_file_existance_and_type,
                        nargs='+',
						help="List of ivar_variants.vcf or ivar_variants.tsv to summarise")
    parser.add_argument('-voc', '--voc_names', type = csv_list,
						 required=True, help="List of Variants of Concern names (e.g. UK, SA, Brazil, Nigeria) "
                        )
    parser.add_argument('-r', '--ref_meta', required=False, default='data/cov_lineage_variants.tsv',
                        type=check_file_existance_and_type,
						help="Path to metadata TSV file containing info on the key mutations")
    parser.add_argument('--signature_snvs_only', required=False, action='store_true',
                        help="Check VCF for only signature/official snvs linked to a VOC")
    parser.add_argument('--stat_filter_snvs', required=False, action='store_true',
                        help="Filter snvs based on statistical significance")
    parser.add_argument('--subplots_orientation',  default="row",
                        required=False, help="List of Variants of Concern names (e.g. UK, SA, Brazil, Nigeria) "
                        )

    args = parser.parse_args()
    output_file_name=""; vcf_df_cleaned=pd.DataFrame()


    VOCmeta_df_full = pd.read_csv(args.ref_meta, sep="\t")
    if args.voc_names == ["all"]:
        vocnames = VOCmeta_df_full["VOC"].unique()
    else:
        vocnames=args.voc_names

    n_heatmaps = len(vocnames) #heatmaps to plot
    n_samples = len(args.input) #samples per plot

    if args.subplots_orientation == "row":
        figsizedim = [1.8*n_heatmaps+(0.5*n_samples),1.7*n_heatmaps] #width and height in inches
        fig, axis = plt.subplots(nrows=1, ncols=n_heatmaps, figsize=figsizedim, dpi=300)
    elif args.subplots_orientation == "column": #column arrangement by default
        figsizedim = [1.8, 1.8*n_heatmaps+(0.5*n_samples)]  # width and height in inches
        fig, axis = plt.subplots(ncols=1, nrows=n_heatmaps, figsize=figsizedim, dpi=300)
    else:
        raise Exception("Subplot orientation values not defined. Permitted values: row, column")

    if not hasattr(axis, '__iter__'): #if number of plots is 1
        axis = [axis]
        fig.set_size_inches(2, 3.5)

    #plots are based on VOCs containing several samples (main VOC plot loop)
    for vocname, subplotaxis in zip(vocnames,axis):

        #filter master metadata based on VOC name parameter
        VOCmeta_df = VOCmeta_df_full[VOCmeta_df_full.VOC == vocname].copy()
        VOCmeta_df.sort_values(by=['Position'], inplace=True)
        VOCmeta_df["NucName+AAName"] = VOCmeta_df["NucName"] + "|" + VOCmeta_df["AAName"]
        nVOCSNVs = VOCmeta_df.shape[0]

        for sample_path in args.input:
            inputtype = get_input_type(sample_path)
            input_file_name = os.path.basename(sample_path)
            output_file_name = str(input_file_name.split(".vcf")[0]) + "." + vocname + ".trimmed.vcf"

            if inputtype == "vcf":
                skip_line_n, vcf_source = find_vcf_headerline_source(sample_path)
                print("VCF source {}\nSelected VOC:{}".format(vcf_source, vocname))
                vcf_df = pd.read_csv(sample_path, sep="\t", skiprows=skip_line_n)
            elif inputtype == "tsv":
                vcf_df = convert_tsv2vcf(sample_path)
                vcf_source = "TSV File"
            else:
                raise Exception("Unsupported input type")

            vcf_df.to_csv("tsv2vcf_temp.vcf", sep="\t",index=False)

            if vcf_df.empty:
                raise  Exception("Empty input vcf_df DataFrame. Input parsing failed")


            if VOCmeta_df.empty:
                raise Exception("Selected VOC {}  is not available in metadata file {}".format(vocname,
                                                                                               args.ref_meta))

            if vcf_source == 'iVar':
                VOCmeta_df.loc[VOCmeta_df["Type"] == "Del", "Position"] += 1


            if args.stat_filter_snvs:
                VOCmeta_df = VOCmeta_df[VOCmeta_df["SignatureSNV"] == True]


            if not type(vcf_df.loc[0,"POS"]) == type(VOCmeta_df.loc[VOCmeta_df.index[0],"Position"]):
                raise Exception("vcf_df[\"POS\"] and VOCmeta_df[\"Position\"] types do not match. Check var type conversions")


            # filter 1: by posistion
            voc_snvs_positions = VOCmeta_df["Position"].to_list()
            vcf_selected_idx = vcf_df["POS"].isin(voc_snvs_positions)



            # filter 2: by match to REF and ALT in metadata for SINGLE BASE Sub and MULTIBASE DELETIONS
            for row_idx_vcf in vcf_df.loc[vcf_selected_idx,:].index:

                #would work as positions in meta and vcf_df match Empty DataFrame should not happen due to position discrep
                Ref_Alt_df = VOCmeta_df[VOCmeta_df["Position"] == vcf_df.loc[row_idx_vcf, "POS"]][["Ref","Alt"]]

                if not Ref_Alt_df.empty:
                    Ref,Alt = Ref_Alt_df.values[0]
                else:
                    raise Exception("No matches for position {} in metadata".format(vcf_df.loc[row_idx_vcf, "POS"]))

                if any(vcf_df.loc[row_idx_vcf, ["REF","ALT"]] == [Ref,Alt]) == False:
                    print("WARNING: Position {} REF and ALT allele mismatch with the metadata. {}/{} (VCF) vs {}/{} (META)".format(
                            vcf_df.loc[row_idx_vcf,"POS"],vcf_df.loc[row_idx_vcf, "REF"],
                            vcf_df.loc[row_idx_vcf, "ALT"],Ref,Alt))

                    vcf_selected_idx[row_idx_vcf]=False #change selection index to false bool

            # filter 3: Add here if present the multi-substitutions positions
            VOCmeta_df_multisub_idx = VOCmeta_df.loc[(VOCmeta_df["Type"] == "Sub") & (VOCmeta_df["Length"] > 1),"Position"].index
            multisub_positions = VOCmeta_df.loc[VOCmeta_df_multisub_idx,"Position"].to_list()

            for idx in VOCmeta_df_multisub_idx:
                multisub_pos = VOCmeta_df.loc[idx,"Position"]
                multisub_len = VOCmeta_df.loc[idx,"Length"]
                multisub_pos_success_counter = [False] * multisub_len #counter of successful n-mer substitutions for decision making
                multisub_pos_success_index = []
                multisub_pos_range_list = list(range(multisub_pos,multisub_pos+multisub_len))


                Ref_meta, Alt_meta = VOCmeta_df.loc[idx, ["Ref", "Alt"]].to_list()
                vcf_df_multisub_pos = vcf_df[vcf_df["POS"].isin(multisub_pos_range_list)] #to increase speed and preserve indices

                #do matches for every position candidate against reference metadata to identify valid mutation sets
                for list_idx, pos in enumerate(multisub_pos_range_list):

                    vcf_df_multisub_pos_idx = vcf_df_multisub_pos.loc[vcf_df_multisub_pos["POS"] == pos].index

                    match_df = vcf_df_multisub_pos[ (vcf_df_multisub_pos["REF"].isin([Ref_meta[list_idx]])) &
                                         (vcf_df_multisub_pos["ALT"].isin([Alt_meta[list_idx]])) ]

                    if not match_df.empty:
                        multisub_pos_success_counter[list_idx] = True
                        multisub_pos_success_index = multisub_pos_success_index + match_df.index.to_list()


                if all(multisub_pos_success_counter):
                    #alt_freq_list = vcf_df.loc[multisub_pos_success_index,
                    #                 vcf_df.columns[-1]].str.split(r':').str[7].astype(float).tolist()
                    #print("Maximum alternative allele frequency for multi-substition is {}".format(max(alt_freq_list)))

                    vcf_df.loc[multisub_pos_success_index[0],["REF","ALT"]] = Ref_meta,Alt_meta
                    vcf_selected_idx[multisub_pos_success_index[0]] = True #want to include only the first position of multi-sub




            vcf_df = vcf_df[vcf_selected_idx]
            VOCmetaNotFound = VOCmeta_df[~VOCmeta_df["Position"].isin(vcf_df["POS"])]

            #DEBUG
            print(VOCmetaNotFound[["VOC","Position","NucName"]])
            print(input_file_name)

            print("In sample {}, a total of {} SNVs were not found:\n {}".format(
                                                   input_file_name,
                                                   VOCmetaNotFound.shape[0],
                                                   VOCmetaNotFound[["NucName", "AAName", "Position"]].to_string(index=False)))

            if len(vcf_df.index) == 0:
                warnings.warn("Zero SNVs found in {} for {} VOC SNVs."
                              "Might be an interesting sample or issue with input ... ".format(input_file_name,
                                                                                               vocname))


            if args.stat_filter_snvs and "FILTER" in vcf_df.columns:
                vcf_df = vcf_df[vcf_df["FILTER"] == "PASS"]

            #filter #4: Remove duplicated entries per position
            vcf_df_cleaned = remove_duplicated_vcf_snvs(vcf_df,VOCmeta_df)
            if vcf_df_cleaned.shape[0] == 0:
                vcf_df_cleaned.loc[0,"CHROM"] = "NO MATCHING SNVs"



            # add extra column for to record sample SNV counts for heatmap
            VOCmeta_df[input_file_name] = [0] * nVOCSNVs
            # append ALT_FREQ values for the selected snvs
            VOCmeta_df.loc[VOCmeta_df["Position"].isin(vcf_df_cleaned["POS"]),input_file_name] = vcf_df_cleaned[vcf_df_cleaned.columns[-1]].str.split(r':').str[7].astype(float).tolist()



            #DEBUG
        VOCmeta_df.to_csv("heatmap_data2plot-{}.tsv".format(vocname),sep="\t")

        if output_file_name and not vcf_df_cleaned.empty:
            print("Writing out {} snvs to VCF".format(vcf_df_cleaned.shape[0]))
            vcf_df_cleaned.to_csv(output_file_name,sep="\t",index=False, mode="w")
            print("Trimmed VCF with VOC snvs is written to \"{}\"".format(output_file_name))


        VOCpangolineage = VOCmeta_df["PangoLineage"].unique()[0]
        VOCheatmapper.renderplot(VOCmeta_df,
                                     title='{} variant ({}) SNVs'.format(vocname, VOCpangolineage),
                                     axis=subplotaxis)




    plt.tight_layout()
    input_folder_name = os.path.basename(os.path.dirname(args.input[0]))
    heatmapfilename="heatmap-overall-{}.png".format(input_folder_name)
    plt.savefig(heatmapfilename)
    print("Heatmap rendered as {} at {}".format(heatmapfilename, os.getcwd()))
    print("Done")
