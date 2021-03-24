#!/usr/bin/env python
import pandas as pd
import argparse, warnings
import os, re
import matplotlib.pyplot as plt

import vcfparser.VOCheatmapper as VOCheatmapper
import vcfparser.BAMutilities as BAMutilities

#constants
support_ext = ["txt","tsv","vcf", "bam"]


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

def check_decimal_range(arg):
    try:
        value = float(arg)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err))

    if value < 0 and value > 1:
        message = "Expected 0 >= value <= 1, got value = {}".format(value)
        raise argparse.ArgumentTypeError(message)

    return value

def main():
    parser = argparse.ArgumentParser("VCFparser.py parses iVar (https://andersen-lab.github.io/ivar/html/manualpage.html) "
                                     "TSV or VCF output files and filters out snvs linked to by the VOC\n")
    parser.add_argument('-i', '--input', required=True,
                        type=check_file_existance_and_type,
                        nargs='+',
						help="List of ivar_variants.vcf or ivar_variants.tsv to summarise")
    parser.add_argument('-bam', '--bam_files', required=False,
                        type=check_file_existance_and_type,
                        nargs='+',
                        help="Optionally provide a list of corresponding bam files in THE SAME ORDER as"
                             "files provided for the -i parameter")
    parser.add_argument('-voc', '--voc_names', type = csv_list,
						 required=True, help="List of Variants of Concern names (e.g. UK, SA, Brazil, Nigeria) "
                        )
    parser.add_argument('-r', '--ref_meta', required=False, default=os.path.dirname(__file__)+'/data/cov_lineage_variants.tsv',
                        type=check_file_existance_and_type,
						help="Path to metadata TSV file containing info on the key mutations")
    parser.add_argument('--signature_snvs_only', required=False, action='store_true',
                        help="Check VCF for only signature/official snvs linked to a VOC")
    parser.add_argument('--key_snvs_only', required=False, action='store_true',
                        help="Check VCF for only the key (S-gene associated)  snvs linked to a VOC")
    parser.add_argument('--stat_filter_snvs', required=False, action='store_true',
                        help="Filter snvs based on statistical significance (i.e. QC PASS flag)")
    parser.add_argument('--subplots_mode',  default="oneplotperfile",
                        required=True, help="How to plot multiple plots (onerow, onecolumn, oneplotperfile)"
                        )
    parser.add_argument('--min_snv_freq_threshold', default=0, type=check_decimal_range, metavar="[0-1]",
                        required=False, help="Set minimum SNV frequency threshold to display (default: 0)"
                        )
    parser.add_argument('--annotate', required=False, action='store_true',
                        help="Annotate heatmap with SNV frequency values")



    args = parser.parse_args()
    output_file_name=""; vcf_df_cleaned=pd.DataFrame()
    input_folder_name = os.path.basename(os.path.dirname(args.input[0]))
    MAXnVOCSNVs=0
    axis_list=[]
    bam_vcf_tsv_files_pairs_dict = {}



    if args.bam_files:
        bam_vcf_tsv_files_pairs_dict = dict(zip(args.input,args.bam_files))


    VOCmeta_df_full = pd.read_csv(args.ref_meta, sep="\t")
    if args.voc_names == ["all"]:
        vocnames = VOCmeta_df_full["VOC"].unique()
    else:
        vocnames=args.voc_names

    if len(vocnames) == 0:
        raise Exception("VOC names not speciefied")

    n_heatmaps = len(vocnames) #heatmaps to plot
    n_samples = len(args.input) #samples per plot


    if args.subplots_mode == "onerow":
        figsizedim = [1.8*n_heatmaps+(0.5*n_samples),1.7*n_heatmaps] #width and height in inches
        fig, axis_list = plt.subplots(nrows=1, ncols=n_heatmaps, figsize=figsizedim, dpi=300)
    elif args.subplots_mode == "onecolumn": #column arrangement by default
        figsizedim = [1.8, 1.8*n_heatmaps+(0.5*n_samples)]  # width and height in inches
        fig, axis_list = plt.subplots(ncols=1, nrows=n_heatmaps, figsize=figsizedim, dpi=300)
    elif args.subplots_mode == "oneplotperfile":
        figures_idx_map_dict={}; fig_counter=1
        for vocname in vocnames:
            fig, axis = plt.subplots(ncols=1, nrows=1, dpi=300)
            figures_idx_map_dict[vocname] = fig_counter
            axis_list.append(axis)
            fig_counter=fig_counter+1
    else:
        raise Exception("Subplot {} orientation is not valid. Permitted values: onerow, onecolumn, oneplotperfile".format(args.subplots_mode))



    if not hasattr(axis_list, '__iter__'): #if number of sub-plots is 1
        axis_list = [axis_list]


    axis_dict = dict(zip(vocnames, axis_list))


    #plots are based on VOCs containing several samples (main VOC plot loop)
    for vocname in vocnames:
        read_coverages_2Darray = []

        #filter master metadata based on VOC name parameter
        VOCmeta_df = VOCmeta_df_full[VOCmeta_df_full.VOC == vocname].copy()
        VOCmeta_df.sort_values(by=['Position'], inplace=True)
        VOCmeta_df["NucName+AAName"] = VOCmeta_df["NucName"] + "|" + VOCmeta_df["AAName"]

        if args.signature_snvs_only:
            VOCmeta_df = VOCmeta_df[VOCmeta_df["SignatureSNV"] == True]
        if args.key_snvs_only:
            VOCmeta_df = VOCmeta_df[VOCmeta_df["Key"] == True]

        nVOCSNVs = VOCmeta_df.shape[0]
        MAXnVOCSNVs = max([MAXnVOCSNVs, nVOCSNVs])

        for sample_path in args.input:
            inputtype = get_input_type(sample_path)
            input_file_name = os.path.splitext(os.path.basename(sample_path))[0] #get samplename without file extension
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


            if not type(vcf_df.loc[0,"POS"]) == type(VOCmeta_df.loc[VOCmeta_df.index[0],"Position"]):
                raise Exception("vcf_df[\"POS\"] and VOCmeta_df[\"Position\"] types do not match. Check var type conversions")


            # filter 1: by posistion
            voc_snvs_positions = VOCmeta_df["Position"].to_list()
            vcf_selected_idx = vcf_df["POS"].isin(voc_snvs_positions)

            # filter 2: by match to REF and ALT alleles in substitutions (SUB) metadata
            #           for SINGLE BASE substituions
            for row_idx_vcf in vcf_df.loc[vcf_selected_idx & (vcf_df["TYPE"] == "SNP"),:].index:

                metadata_pos_idx = VOCmeta_df["Position"] == vcf_df.loc[row_idx_vcf, "POS"]

                #would work as positions in meta and vcf_df match Empty DataFrame should not happen due to position discrep
                Ref_Alt_df = VOCmeta_df[ metadata_pos_idx ][["Ref","Alt"]]


                if not Ref_Alt_df.empty:
                    Ref,Alt = Ref_Alt_df.values[0]
                else:
                    raise Exception("No matches for position {} in metadata".format(vcf_df.loc[row_idx_vcf, "POS"]))


                if not all(vcf_df.loc[row_idx_vcf, ["REF","ALT"]] == [Ref,Alt]):
                    print("WARNING: Position {} REF and ALT allele mismatch with the metadata. {}/{} (VCF) vs {}/{} (META)".format(
                            vcf_df.loc[row_idx_vcf,"POS"],vcf_df.loc[row_idx_vcf, "REF"],
                            vcf_df.loc[row_idx_vcf, "ALT"],Ref,Alt))

                    vcf_selected_idx[row_idx_vcf]=False #change selection index to false bool

            #filter Deletions separately just by position and deletion length
            #makes compatible with iVar and Virontus pipeline formats
            for row_idx_vcf in vcf_df.loc[vcf_selected_idx & (vcf_df["TYPE"] == "DEL"), :].index:
                metadata_pos_idx = VOCmeta_df[VOCmeta_df["Position"] == vcf_df.loc[row_idx_vcf, "POS"]].index

                if len(metadata_pos_idx) >= 2:
                    raise Exception("Metadata should have unique deletion definitions."
                                    "Offending value {} for VOC {}"
                                    .format(VOCmeta_df.loc[metadata_pos_idx,'NucName+AAName'].to_list()[0],
                                            vocname))
                if (int(VOCmeta_df.loc[metadata_pos_idx,"Length"])+1 != len(vcf_df.loc[row_idx_vcf,"REF"])):
                    vcf_selected_idx[row_idx_vcf] = False
                    raise Warning("Deletion length at position {} did not match. Skipping this position ...".format(int(VOCmeta_df[metadata_pos_idx]["Position"])))


            # filter 3: Add here if present the multi-substitutions positions
            VOCmeta_df_multisub_idx = VOCmeta_df.loc[(VOCmeta_df["Type"] == "Sub") & (VOCmeta_df["Length"] > 1),"Position"].index
            #multisub_positions = VOCmeta_df.loc[VOCmeta_df_multisub_idx,"Position"].to_list()

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




            vcf_df_temp = vcf_df[vcf_selected_idx].copy()
            VOCmetaNotFound = VOCmeta_df[~VOCmeta_df["Position"].isin(vcf_df_temp["POS"])]

            #DEBUG
            #print(VOCmetaNotFound[["VOC","Position","NucName"]])
            #print(input_file_name)

            print("In sample {}, a total of {} SNVs were not found:\n {}".format(
                                                   input_file_name,
                                                   VOCmetaNotFound.shape[0],
                                                   VOCmetaNotFound[["NucName", "AAName", "Position"]].to_string(index=False)))

            if len(vcf_df_temp.index) == 0:
                warnings.warn("Zero SNVs found in sample {} for {} VOC SNVs."
                              "Might be an interesting sample or issue with input ... ".format(input_file_name,
                                                                                               vocname))


            if args.stat_filter_snvs and "FILTER" in vcf_df_temp.columns:
                vcf_df_temp = vcf_df_temp[vcf_df_temp["FILTER"] == "PASS"]

            #filter #4: Remove duplicated entries per position
            vcf_df_cleaned = remove_duplicated_vcf_snvs(vcf_df_temp,VOCmeta_df)
            if vcf_df_cleaned.shape[0] == 0:
                vcf_df_cleaned.loc[0,"CHROM"] = "NO MATCHING SNVs"



            # ASSIGMENT OF FREQ VALUES
            # add extra column for to record sample SNV ALT_FREQ for heatmap
            VOCmeta_df[input_file_name] = [0] * nVOCSNVs

            vcf_df_cleaned["ALT_FREQ"] = vcf_df_cleaned[vcf_df_cleaned.columns[-1]].str.split(r':').str[7].astype(float).tolist()
            VOCmeta_df.input_file_name = 0

            if all(vcf_df_cleaned["POS"].isna()):
                warnings.warn("NO MATCHING SNVs were found for VOC {}. Skipping this VOC ... ".format(vocname))
                continue



            for idx in vcf_df_cleaned.index:
                if vcf_df_cleaned.loc[idx,"TYPE"]=="SNP":
                    match_query = 'Position == '+ str(vcf_df_cleaned.loc[idx,"POS"])+\
                                  ' & Ref == \"'+vcf_df_cleaned.loc[idx,"REF"] +'\"'+\
                                  ' & Alt == \"'+vcf_df_cleaned.loc[idx,"ALT"]+'\"''                                        '
                else:
                    match_query = 'Position == ' + str(vcf_df_cleaned.loc[idx, "POS"])

                VOCmeta_sel_index = VOCmeta_df.query(match_query).index
                if not VOCmeta_sel_index.empty:
                    VOCmeta_df.loc[VOCmeta_sel_index,input_file_name] = vcf_df_cleaned.loc[idx,"ALT_FREQ"]




            #filter based on set threshold if set by the user
            if args.min_snv_freq_threshold:
                VOCmeta_df.loc[VOCmeta_df[input_file_name] <= args.min_snv_freq_threshold, input_file_name] = 0

            #check positions for coverage
            if bam_vcf_tsv_files_pairs_dict:
                read_coverages_2Darray.append(BAMutilities.get_ref_coverage_array(bam_vcf_tsv_files_pairs_dict[sample_path],
                                                             query_positions=VOCmeta_df["Position"]))




        #DEBUG
        VOCmeta_df.sort_values("Position",inplace=True)
        VOCmeta_df.to_csv("heatmap_data2plot-{}.tsv".format(vocname),sep="\t")
        print("INFO: Data to plot written to heatmap_data2plot-{}.tsv".format(vocname))

        if output_file_name and not vcf_df_cleaned.empty:
            print("INFO: Writing out {} snvs to VCF".format(vcf_df_cleaned.shape[0]))
            vcf_df_cleaned.to_csv(output_file_name,sep="\t",index=False, mode="w")
            print("INFO: Trimmed VCF with VOC snvs is written to \"{}\"".format(output_file_name))


        #render plot
        VOCpangolineage = VOCmeta_df["PangoLineage"].unique()[0]
        if args.subplots_mode == "oneplotperfile":
            if nVOCSNVs < 6:
                ysizein = nVOCSNVs*0.135+2
            else:
                ysizein = nVOCSNVs * 0.135

            fig = plt.figure(figures_idx_map_dict[vocname]) #make figure object active for rendering
            fig.set_size_inches(1.8+0.1*n_samples,ysizein)
            VOCheatmapper.renderplot(VOCmeta_df,
                                     title='{} variant ({}) SNVs'.format(vocname, VOCpangolineage),
                                     axis=axis_dict[vocname],
                                     is_plot_annotate=args.annotate,
                                     read_coverages_2Darray=read_coverages_2Darray)
            heatmapfilename = "heatmap-overall-{}-{}.png".format(input_folder_name, vocname)
            plt.tight_layout()
            plt.savefig(heatmapfilename)
            plt.close()
            print("INFO: Heatmap rendered as {} at {}".format(heatmapfilename, os.getcwd()))
        else:
            VOCheatmapper.renderplot(VOCmeta_df,
                                     title='{} variant ({}) SNVs'.format(vocname, VOCpangolineage),
                                     axis=axis_dict[vocname],
                                     is_plot_annotate=args.annotate)







    #render subplots if this feature is selected
    if plt.get_fignums(): #if any figures left open, then save them to file
        fig.set_size_inches(fig.get_size_inches()[0], 0.10 * MAXnVOCSNVs)  # width and height
        plt.tight_layout()
        heatmapfilename="heatmap-overall-{}.png".format(input_folder_name)
        plt.savefig(heatmapfilename)
        plt.close()
        print("INFO: Heatmap rendered as {} at {}".format(heatmapfilename, os.getcwd()))

    print("Done")

if __name__ == '__main__':
    main()
#TODO: bold S-gene linked SNVs in heatmap
#TODO: for non-called snvs (i.e. with no-frequency) check for coverage and if coverage is absent annotate as NoCov or NC
#TODO: Add colour range (0,0.1) while non-called SNVs will be white - Done
#TODO: remove tsv extension from the sample names - Done
#TODO: add optional frequencies text annotation key
#TODO: make y-axis more squeezed (less space between snv names) - Done

#TODO: sort the heatmap by position as data frame manipulation randomizes snvs (y-axis)
#TODO: make color bar legend independent of y-axis length
#TODO: input file with input tsv and bam file pairs and sample name
#TODO: make deletions detections independent of the  actual deleted sequence. Only check for position and length only. Allow accept both viralrecon and iVar