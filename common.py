import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
from collections import defaultdict
import pickle


# function to read in and format the fusion master
# new fusion master is /nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/fusion_master_20240415.txt
# old is /nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/fusion_master.txt
def read_fm(fm_path):
    # if new fusion master...
    if fm_path=='/nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/fusion_master_20240415.txt':
        fm = pd.read_csv(fm_path, sep='\t', na_filter=False)
        # filter to just b-all
        fm["b-all_count"] = fm["b-all_count"].replace("", 0)
        fm["b-all_count"] = fm["b-all_count"].astype(int)
        fm = fm.loc[fm["b-all_count"]>1, :]
        # create gene1 and gene2 columns
        fm['gene1'] = fm.gene
        # convert the gene partner column into list
        fm["gene2"] = fm["b-all_fusion_partners"].apply(lambda x: x.split('; '))
        # explode the gene partner column into multiple rows
        fm = fm.explode('gene2')
        # single gene or has gene partner?
        fm['single']='N'
        fm.loc[fm['gene2'] == '','single'] = 'Y'
        fm = fm[['gene1', 'gene2', 'single']]
        # create different permutations
        fm.loc[:, 'fusion_std12'] = None
        fm.loc[:, 'fusion_std21'] = None
        fm = fm.reset_index(drop=True)
        fm.loc[fm['single']=='N', "fusion_std12"] = fm['gene1'] + "::" + fm['gene2']
        fm.loc[fm['single']=='N', "fusion_std21"] = fm['gene2'] + "::" + fm['gene1']
        # drop dupes 
        fm['dupe_id'] = fm[['gene1', 'gene2']].values.tolist()
        fm.dupe_id = fm.dupe_id.sort_values().apply(lambda x: sorted(x))
        fm = fm[~fm.dupe_id.duplicated()]
        fm = fm.drop("dupe_id",axis=1)
        fm = fm.reset_index(drop=True)
    # if old fusion master...
    else:
        # read in gene fusion master
        fm = pd.read_csv(fm_path, sep='\t')
        # reformat the master list
        fm[['gene1', 'gene2']] = fm['fusion'].str.split('-', n=1, expand=True)
        fm['single'] = 'N'
        fm.loc[fm['gene2'].isnull(), 'single'] = 'Y'
        fm.loc[:, 'fusion_std12'] = None
        fm.loc[:, 'fusion_std21'] = None
        fm.loc[fm['single']=='N', "fusion_std12"] = fm['gene1'] + "::" + fm['gene2']
        fm.loc[fm['single']=='N', "fusion_std21"] = fm['gene2'] + "::" + fm['gene1']
    return fm

def std_lookup(val, std_lookup_path):
    std = pd.read_csv(std_lookup_path, sep='\t', index_col=0)
    std = std.to_dict('dict')['std_name']
    if val in std:
        return std[val]
    else:
        print(f"{val} is not in {std_lookup_path}")
        
# function to read and format the metadata file
def read_md(md_path, fusion_std_path, gu_std_fusion_path, gu_pax5_status_std_path):
    print("\nReading in metadata...")
    # read in metadata
    # read in metadata file
    metadata = pd.read_csv(md_path, sep='\t')
    # drop any rows that are null
    metadata = metadata[metadata['seq_id'].notna()]
    # read in dictionary for standardizing fusion names...this list will need to be updated with new fusions
    # since a fusion partner gene can be promiscuous ie CRLF2, the std_name is read in as a list
    # fusion_std = pd.read_csv(fusion_std_path, sep='\t', index_col=0, converters={"std_name": eval})
    fusion_std = pd.read_csv(fusion_std_path, sep='\t', index_col=0)
    gu_std = pd.read_csv(gu_std_fusion_path, sep='\t', index_col=0)
    pax5_std = pd.read_csv(gu_pax5_status_std_path, sep='\t', index_col=0)
    ges_fusion_yn = fusion_std.to_dict('dict')['fusion']
    fusion_std2 = fusion_std.to_dict('dict')['std_name2']
    fusion_std = fusion_std.to_dict('dict')['std_name']
    gu_fusion_yn = gu_std.to_dict('dict')['fusion']
    gu_std = gu_std.to_dict('dict')['std_name']
    gu_p5_fusion_yn = pax5_std.to_dict('dict')['fusion']
    pax5_std = pax5_std.to_dict('dict')['std_name']
    metadata['genomic_subtype_std'] = metadata.genomic_subtype.map(fusion_std)
    metadata['ges_fusion_yn'] = metadata.genomic_subtype.map(ges_fusion_yn)
    metadata['genomic_subtype_std2'] = metadata.genomic_subtype.map(fusion_std2)
    metadata['gu_fusion_std'] = metadata.gu_fusion.map(gu_std)
    metadata['gu_fusion_yn'] = metadata.gu_fusion.map(gu_fusion_yn)
    metadata['gu_pax5_status_std'] = metadata.gu_pax5_status.map(pax5_std)
    metadata['gu_pax5_fusion_yn'] = metadata.gu_pax5_status.map(gu_p5_fusion_yn)
    # metadata.loc[metadata.gu_fusion_std.isna(), 'gu_fusion_std'] = '' #?? should we do unknown or something else
    # metadata.loc[metadata.gu_pax5_status_std.isna(), 'gu_pax5_status_std'] = ''#?? should we do unknown or something else
    # create new row for seq_ids with promiscuous partner genes as those partner genes will have more than one genomic_subtype_std
    # ie Ph-like CRLF2r --> ['IGH::CRLF2', 'P2RY8::CRLF2']
    # metadata = metadata.explode('genomic_subtype_std')
    # create column that summarizes all the subtype info
    # metadata.loc[:, 'genomic_subtype_std3'] = metadata.genomic_subtype_std
    # metadata.loc[metadata.genomic_subtype_std3 == 'unknown', 'genomic_subtype_std2'] = None
    # metadata.loc[metadata.genomic_subtype_std3 == 'none', 'genomic_subtype_std2'] = None
    # metadata.loc[:, 'gu_fusion_std2'] = metadata.gu_fusion_std
    # metadata.loc[metadata.gu_fusion_std2 == 'unknown', 'gu_fusion_std2'] = None
    # metadata.loc[metadata.gu_fusion_std2 == 'none', 'gu_fusion_std2'] = None
    # metadata.loc[:, 'gu_pax5_status_std2'] = metadata.gu_pax5_status_std
    # metadata.loc[metadata.gu_pax5_status_std2 == 'unknown', 'gu_pax5_status_std2'] = None
    # metadata.loc[metadata.gu_pax5_status_std2 == 'none', 'gu_pax5_status_std2'] = None
    
    # metadata.loc[:, 'genomic_subtype_agg_std'] = ''
    # for index, row in metadata.iterrows():
    #     # make unique list
    #     metadata.at[index, 'genomic_subtype_agg_std'] = list(dict.fromkeys([row.genomic_subtype_std3, row.gu_fusion_std2, row.gu_pax5_status_std2]))
    #     # get rid of None
    #     metadata.at[index, 'genomic_subtype_agg_std'] = [x for x in metadata.at[index, 'genomic_subtype_agg_std'] if x is not None]
    #     # convert to string
    #     metadata.at[index, 'genomic_subtype_agg_std'] = '|'.join(metadata.at[index, 'genomic_subtype_agg_std'])
    #     # convert to list and dedupe
    #     metadata.at[index, 'genomic_subtype_agg_std'] = metadata.at[index, 'genomic_subtype_agg_std'].split("|")
    #     metadata.at[index, 'genomic_subtype_agg_std'] = list(dict.fromkeys(metadata.at[index, 'genomic_subtype_agg_std']))
    #     # convert back to string
    #     metadata.at[index, 'genomic_subtype_agg_std'] = '|'.join(metadata.at[index, 'genomic_subtype_agg_std'])
        # list(dict.fromkeys(mylist))
        # print(row.genomic_subtype_agg_std)
    # metadata.loc[:, 'genomic_subtype_agg_std'] = list(metadata.genomic_subtype_std2, metadata.gu_fusion_std2, metadata.gu_pax5_status_std2)
    
    # metadata = metadata.drop(['genomic_subtype_std3', 'gu_fusion_std2', 'gu_pax5_status_std2'], axis = 1)
    
    # warnings
    for std_type in [('genomic_subtype', 'genomic_subtype_std'), 
                     ('gu_fusion', 'gu_fusion_std'), 
                     ('gu_pax5_status', 'gu_pax5_status_std')]:
        if metadata[std_type[1]].isna().any():
            print(f"WARNING: These are the {std_type[0]} values that are NOT standardized!:\n")
            print(metadata[std_type[0]][metadata[std_type[1]].isna()].value_counts())
        else:
            print(f"Everything is standardized for {std_type[1]}!")
    
    # check discrepancies...
    diff = metadata.loc[(metadata.ges_fusion_yn == 'Y') & (metadata.gu_fusion_std == 'none'), 
                        ['seq_id', 'genomic_subtype_std', 'gu_fusion_std']]
    if diff.shape[0] != 0:
        print("WARNING: There are some metadata discrepancies with fusions reported from Gu:")
        print(diff)
    diff = metadata.loc[(metadata.ges_fusion_yn == 'N') & (metadata.gu_fusion_yn == 'Y'),
                        ['seq_id', 'genomic_subtype_std', 'gu_fusion_std']]
    if diff.shape[0] != 0:
        print("WARNING: There are some metadata discrepancies with fusions reported from Gu:")
        print(diff)
    diff = metadata.loc[(metadata.ges_fusion_yn == 'Y') & (metadata.gu_pax5_status_std == 'none'),
                        ['seq_id', 'genomic_subtype_std', 'gu_pax5_status_std']]
    if diff.shape[0] != 0:
        print("WARNING: There are some metadata discrepancies with fusions reported from Gu:")
        print(diff)
    diff = metadata.loc[(metadata.ges_fusion_yn == 'N') & (metadata.gu_pax5_fusion_yn == 'Y'),
                        ['seq_id', 'genomic_subtype_std', 'gu_pax5_status_std']]
    if diff.shape[0] != 0:
        print("WARNING: There are some metadata discrepancies with fusions reported from Gu:")
        print(diff)
    diff = metadata.loc[(metadata.gu_fusion_yn == 'N') & (metadata.gu_pax5_fusion_yn == 'Y'),
                        ['seq_id', 'gu_fusion_std', 'gu_pax5_status_std']]
    if diff.shape[0] != 0:
        print("WARNING: There are some metadata discrepancies with fusions reported from Gu:")
        print(diff)

    for std_type in [('genomic_subtype', 'genomic_subtype_std', 'ges_fusion_yn'), 
                     ('gu_fusion', 'gu_fusion_std', 'gu_fusion_yn'), 
                     ('gu_pax5_status', 'gu_pax5_status_std', 'gu_pax5_fusion_yn')]:
        print(f"For {std_type[0]} these are the standardizations in your data:\n")
        print(metadata[[std_type[0], std_type[1], std_type[2]]].value_counts(dropna=False).reset_index().sort_values(by=['count', std_type[0]], ascending=False))
    return metadata

# function to compare the fusions in fusion_std and fusion_master
def compare_fstd_fm(fusion_std_path, fm):
    # read in fusion_std
    # fusion_std = pd.read_csv(fusion_std_path, sep='\t', index_col=0, converters={"std_name": eval})
    # fusion_std = fusion_std.explode('std_name')
    fusion_std = pd.read_csv(fusion_std_path, sep='\t', index_col=0)    
    print("\nComparing fusion_std and fusion_master...")
    fs = []
    for i1 in set(fusion_std.std_name):
        if '|' in i1:
            for i2 in i1.split('|'):
                fs.append(i2)
        else:
            fs.append(i1)
    fm = list(fm.fusion_std21) + list(fm.fusion_std12) + list(fm['gene1'][fm.single == "Y"])
    diff = set(fs) - set(fm)
    if len(diff) != 0:
        print("WARNING: Fusions in fusion_std NOT in fusion_master! (maybe these are T-ALL/AML fusions):\n")
        print(diff)
    
# function to compare the fusions in metadata and fusion_master
def compare_md_fm(metadata, fm):
    print("\nComparing fusions in metadata and fusion_master...")
    fm = list(fm.fusion_std21) + list(fm.fusion_std12) + list(fm['gene1'][fm.single == "Y"])
    for std_type in [('genomic_subtype_std', 'ges_fusion_yn'), ('gu_fusion_std', 'gu_fusion_yn'), ('gu_pax5_status_std', 'gu_pax5_fusion_yn')]:
        mfus = list(np.unique(metadata.loc[(metadata[std_type[1]] == 'Y'), std_type[0]]))
        mfus2 = []
        for i1 in set(mfus):
            if '|' in i1:
                for i2 in i1.split('|'):
                    mfus2.append(i2)
            else:
                mfus2.append(i1)
        diff = set(mfus2) - set(fm)
        if len(diff) > 0:
            print(f"WARNING: There are fusions from {std_type[0]} in your metadata NOT in fusion_master!\n{diff}")
            print("If the partner genes above are not in the peds.leuk.bed.XXXXXX file...you can ignore them. Also, these may be T-ALL/AML fusions, which you can also ignore.")
    # assert(len(diff) == 0), "WARNING: There are fusions in your metadata NOT in fusion_master!\n"+str(diff)
    
# function to read in results of given algorithm categorize results --> TPs, FPs, TNs, FNs
# and output the algorithm's performance metrics
# fm is the fusion master as read in by read_fm
def def_results(df, metadata, fm, name, fd_output_path, metadata_fusion, fusion_std_path, gu_fusion_std_path, gu_pax5_status_std_path, ls=['B-ALL'], wlist_path=''):
    if wlist_path != 'none':
        wlist_tf = True
        # read in the whitelist file
        wlist = []
        with open(wlist_path) as f:
            next(f)
            for line in f:
                wlist.append(line.strip())
    else:
        wlist_tf = False
    print("\nStarting evaluation of F1 stats...")
    df = df.merge(metadata[['seq_id', 'lineage', 'genomic_subtype', 'genomic_subtype_std', 'genomic_subtype_std2', 'ges_fusion_yn', 
                            'gu_fusion', 'gu_fusion_std', 'gu_fusion_yn',
                            'gu_pax5_status', 'gu_pax5_status_std', 'gu_pax5_fusion_yn']], on='seq_id')


    num_seqids = len(np.unique(df.seq_id))


    # filter for just desired lineages (ls)
    if (len(ls) != 0) and (metadata_fusion == True):
        print(f"Filtering for lineages {ls}...")
        l0 = len(np.unique(df.seq_id))
        df = df[df.lineage.isin(ls)]
        l1 = len(np.unique(df.seq_id))
        print(f"There were {l0-l1} seq_ids that were removed due to lineage filtering (ie B-ALL).")
    else:
        l0 = len(np.unique(df.seq_id))
        l1 = l0

    # filter out any unknowns...
    if metadata_fusion == True:
        print("Filtering out unknowns...")
        u0 = len(np.unique(df.seq_id))
        df = df.loc[df.genomic_subtype_std != 'unknown', :] #??
        u1 = len(np.unique(df.seq_id))
        print(f"There were {u0-u1} seq_ids that had unknown genomic_std_subtypes and were removed from consideration.")
    else:
        u0 = len(np.unique(df.seq_id))
        u1 = u0

    # compile master list of fusions
    fm_list = list(fm['fusion_std12'][fm['single'] == 'N']) + list(fm['fusion_std21'][fm['single'] == 'N'])
    # there are no single fusions in the B-ALL set....so likely can get rid of this??
    fm_list_single = list(fm['gene1'][fm['single'] == 'Y'])
    # fm_list_single_strings = '|'.join(fm_list_single)

    # create column to denote relevant fusions detected
    print("Marking relevant fusions...")
    df.loc[:, 'fusion_detected_rel'] = 'N'
    df.loc[df.fusion_detected_std12.isin(fm_list) | 
           df.fusion_detected_std21.isin(fm_list) | 
           df.gene1.isin(fm_list_single) | 
           df.gene2.isin(fm_list_single), 'fusion_detected_rel'] = 'Y'

    # create columns to denote relevant fusions in metadata
    df.loc[:, 'genomic_subtype_std_rel'] = ''
    df.loc[:, 'gu_fusion_std_rel'] = ''
    df.loc[:, 'gu_pax5_status_std_rel'] = ''
    # create list of seq_ids with a relevant fusion across any metadata
    all_tp_seqdids = []
    all_tp_seqdids_std_lkup = {}
    for index, row in df.iterrows():
        for std_type in ['genomic_subtype_std', 'gu_fusion_std', 'gu_pax5_status_std']:
            yn_str = []
            for f in row[std_type].split("|"):
                if "::" in f:
                    if f in fm_list:
                        yn_str.append("Y")
                        all_tp_seqdids.append(row.seq_id)
                        if row.seq_id not in all_tp_seqdids_std_lkup:
                            all_tp_seqdids_std_lkup[row.seq_id] = {f}
                        else:
                            all_tp_seqdids_std_lkup[row.seq_id].add(f)
                    elif (f.split("::")[0] in fm_list_single) or (f.split("::")[1] in  fm_list_single):
                        yn_str.append("Y")
                        all_tp_seqdids.append(row.seq_id)
                        if row.seq_id not in all_tp_seqdids_std_lkup:
                            all_tp_seqdids_std_lkup[row.seq_id] = {f}
                        else:
                            all_tp_seqdids_std_lkup[row.seq_id].add(f)
                    else:
                        yn_str.append("N")
                elif f in fm_list_single:
                    yn_str.append("Y")
                    all_tp_seqdids.append(row.seq_id)
                    if row.seq_id not in all_tp_seqdids_std_lkup:
                            all_tp_seqdids_std_lkup[row.seq_id] = {f}
                    else:
                        all_tp_seqdids_std_lkup[row.seq_id].add(f)
                else:
                    yn_str.append("N")
            df.loc[index, std_type+'_rel'] = "|".join(yn_str)
    all_tp_seqdids = list(np.unique(all_tp_seqdids))
    
    all_tp_seqdids_lkup = {}
    cols = ['genomic_subtype', 'gu_fusion', 'gu_pax5_status']
    subset = df.loc[df.seq_id.isin(all_tp_seqdids), ['seq_id'] + cols].drop_duplicates()
    for col in cols:
        for index, row in subset.iterrows():
            if row.seq_id not in all_tp_seqdids_lkup:
                all_tp_seqdids_lkup[row.seq_id] = [row[col]]
            else:
                all_tp_seqdids_lkup[row.seq_id].append(row[col])
    all_tp_seqdids_lkup

    print("evaluating true positives...")
    df['TP_fusion'] = 0
    if metadata_fusion == True:
        # standard fusion gene targets ie gene1::gene2
        df.loc[(df.fusion_detected_rel == 'Y') & 
               ((df.fusion_detected_std12 == df.genomic_subtype_std) |
                (df.fusion_detected_std21 == df.genomic_subtype_std) |
                (df.fusion_detected_std12 == df.gu_fusion_std) |
                (df.fusion_detected_std21 == df.gu_fusion_std) |
                (df.fusion_detected_std12 == df.gu_pax5_status_std) |
                (df.fusion_detected_std21 == df.gu_pax5_status_std)), 'TP_fusion'] = 1

        # for index, row in df.iterrows():
        #     if row.fusion_detected_rel == 'Y':
        #         for std_type in ['genomic_subtype_std', 'gu_fusion_std', 'gu_pax5_status_std']:
        #             for f in row[std_type].split("|"):
        #                 if "::" in f:
        #                     if (row.fusion_detected_std12 == f) | (row.fusion_detected_std21 == f):
        #                         df.loc[index, 'TP_fusion'] = 1
        #                     elif ((f.split("::")[0] in fm_list_single) & (row.gene1 == f.split("::")[0])) | ((f.split("::")[1] in fm_list_single) & (row.gene1 == f.split("::")[1])) | ((f.split("::")[0] in fm_list_single) & (row.gene2 == f.split("::")[0])) | ((f.split("::")[1] in fm_list_single) & (row.gene2 == f.split("::")[1])):
        #                         df.loc[index, 'TP_fusion'] = 1
        #                 elif (f in fm_list_single) & ((row.gene1 == f) | (row.gene2 == f)):
        #                     df.loc[index, 'TP_fusion'] = 1

        for index, row in df.iterrows():
            if row.fusion_detected_rel == 'Y':
                for std_type in ['genomic_subtype_std', 'gu_fusion_std', 'gu_pax5_status_std']:
                    for f in row[std_type].split("|"):
                        if "::" in f:
                            if (row.fusion_detected_std12 == f) | (row.fusion_detected_std21 == f):
                                df.loc[index, 'TP_fusion'] = 1
                            elif ((f.split("::")[0] in fm_list_single) & (row.gene1 == f.split("::")[0])) | ((f.split("::")[1] in fm_list_single) & (row.gene1 == f.split("::")[1])) | ((f.split("::")[0] in fm_list_single) & (row.gene2 == f.split("::")[0])) | ((f.split("::")[1] in fm_list_single) & (row.gene2 == f.split("::")[1])):
                                df.loc[index, 'TP_fusion'] = 1
                        elif (f in fm_list_single) & ((row.gene1 == f) | (row.gene2 == f)):
                            df.loc[index, 'TP_fusion'] = 1 
        df.loc[:, 'TP_seq_id'] = 0
        for seq_id in all_tp_seqdids_lkup:
            # list of fusions in unstandard form
            fs = all_tp_seqdids_lkup[seq_id]
            # number of fusions to detect
            numf = len(fs)
            # establish counter of fusions detected
            numd = 0
            # subset to the fusions detected for given seq_Id
            df_sub = df.loc[df.seq_id == seq_id]
            # type_lkup = {0:'gs', 1:'gu', 2:'gu_pax5'}
            typ = 0
            # iterate through list of fusions
            for f in fs:
                if pd.isnull(f):
                    typ += 1
                    numf = numf - 1
                    continue
                else:
                    if typ == 0:
                        fstds = std_lookup(f, fusion_std_path).split("|")
                    elif typ == 1:
                        fstds = std_lookup(f, gu_fusion_std_path).split("|")
                    elif typ == 2:
                        fstds = std_lookup(f, gu_pax5_status_std_path).split("|")
                    typ += 1
                    # if there is only 1 target standard fusions ie ETV-RUNX1 --> ETV6::RUNX1...check if that was detected
                    if len(fstds) == 1:
                        fstd = fstds[0]
                        if fstd not in fm_list:
                            numf = numf - 1
                        else:
                            if (fstd in df_sub.fusion_detected_std12.tolist()) | (fstd in df_sub.fusion_detected_std21.tolist()):
                                numd += 1
                    # if there is > 1 standard fusions ie PAX5r --> PAX5::AUTS2|PAX5::DACH1......check if any was detected
                    elif len(fstds) > 1:
                        # some unstandard values have ; in them indicating multiple fusions to find ie all of them must be found
                        if (';' in f) & (len(f.split(';')) == len(fstds)):
                            # increment the number of fusions to detect based on number of separators
                            fadd = f.count(';')
                            numf = numf + fadd
                            for fstd in fstds:
                                if fstd not in fm_list:
                                    numf = numf - 1
                                else:
                                    if (fstd in df_sub.fusion_detected_std12.tolist()) | (fstd in df_sub.fusion_detected_std21.tolist()):
                                        numd += 1
                        # others represented possible rearrangement in which case any of them found would suffice
                        else:
                            flag = 0
                            # get number of fusions to detect
                            numf2 = len(fstds)
                            # first remove any fusions not in fusion master
                            for fstd in fstds:
                                if fstd not in fm_list:
                                    fstds.remove(fstd)
                                    numf2 = numf2 - 1
                            # if none of the fusions are in fm_list then we can decrease the number of things to detect
                            if numf2 ==0:
                                numf = numf - 1
                                continue
                            for fstd in fstds:
                                if (fstd in df_sub.fusion_detected_std12.tolist()) | (fstd in df_sub.fusion_detected_std21.tolist()):
                                    flag = 1
                            if flag == 1:
                                numd += 1
                    
            if numf == numd:
                df.loc[df.seq_id == seq_id, 'TP_seq_id'] = 1
            
    
                    
                        # else:
                        #     yn_str.append("N")
    
    # for std_type in ['genomic_subtype_std', 'gu_fusion_std', 'gu_pax5_status_std']:
    #     # need another condition for single genes
    #     df.loc[(df.fusion_detected_rel == 'Y') &
    #            (df[std_type].isin(fm_list_single)) &
    #            ((df.gene1 == df[std_type]) |
    #             (df.gene2 == df[std_type])), 'TP'] = 1
    #     # need another condition for genomic_subtype_std that has multiple targets
    #     for index, row in df.iterrows():
    #         if (row.fusion_detected_rel == 'Y') & ("|" in row[std_type]):
    #             for f in row[std_type].split("|"):
    #                 # double partner genes
    #                 if (row.fusion_detected_std12 == f) | (row.fusion_detected_std21 == f):
    #                     df.loc[index, 'TP'] = 1
    #                 # single genes
    #                 if (f in fm_list_single) & ((row.gene1 == f) | (row.gene2 == f)):
    #                     df.loc[index, 'TP'] = 1

        tps = np.unique(df.seq_id[(df.TP_seq_id == 1)])
        tps = list(set(tps))
        num_tp = len(tps)

    # false negatives
    print("evaluating false negatives...")
    df['FN'] = 0
    fns = []
    if metadata_fusion == True:
        for index, row in df.iterrows():
            for std_ty_rel in ['genomic_subtype_std_rel', 'gu_fusion_std_rel', 'gu_pax5_status_std_rel']:
                for rel in row[std_ty_rel].split("|"):
                    if rel == "Y":
                        fns.append(row.seq_id)
        fns = np.unique(fns)
        # which samples have a genomic subtype but ARE not marked as true positive...
        fns = list(set(fns).difference(tps))
        num_fn = len(fns)
        # mark the false negative
        df.loc[df['seq_id'].isin(fns), "FN"] = 1
    
    # false positives
    print("evaluating false positives...")
    df['FP_fusion'] = 0
    df['FP_seq_id'] = 0
    if metadata_fusion == True:
        df.loc[(df.fusion_detected_rel == 'Y') & (df.TP_fusion != 1), 'FP_fusion'] = 1
        for index, row in df.iterrows():
            # marker for whether the sample has a B-ALL relevant subtype in our metadata
            rel_final = "N"
            for std_ty_rel in ['genomic_subtype_std_rel', 'gu_fusion_std_rel', 'gu_pax5_status_std_rel']:
                for rel in row[std_ty_rel].split("|"):
                    if rel == "Y":
                        rel_final = "Y"
            if wlist_tf == False:
                if (rel_final == 'N') & (row.fusion_detected_rel == "Y"):
                    df.loc[df.seq_id == row.seq_id, 'FP_seq_id'] = 1
            elif wlist_tf == True:
                if (rel_final == 'N') & (row.fusion_detected_rel == "Y"):
                    if (row.gene1+'::'+row.gene2 not in wlist) & (row.gene2+'::'+row.gene1 not in wlist):
                        df.loc[df.seq_id == row.seq_id, 'FP_seq_id'] = 1
                    # turn FP_seq_id col back to 0 for those detected fusions that are in whitelist
                    elif ((row.gene1+'::'+row.gene2 in wlist) | (row.gene2+'::'+row.gene1 in wlist)) & (row.FP_fusion == 1):
                        df.loc[(df.seq_id == row.seq_id) & (df.gene1 == row.gene1) & (df.gene2 == row.gene2), 'FP_seq_id'] = 0
        fps = np.unique(df.seq_id[(df.FP_seq_id == 1)])
        fps = list(set(fps))
        num_fp = len(fps)
    
    # true negatives
    print("evaluating true negatives...")
    df['TN'] = 0
    if metadata_fusion == True:
        # has no fusion genomic_subtype_std and none of the fusion detected are relevant ie it is not a false positive
        # df.loc[(df.ges_fusion_yn == 'N') & (~df.seq_id.isin(fps)), 'TN'] = 1
        for index, row in df.iterrows():
            rel_final = "N"
            for std_ty_rel in ['genomic_subtype_std_rel', 'gu_fusion_std_rel', 'gu_pax5_status_std_rel']:
                for rel in row[std_ty_rel].split("|"):
                    if rel == "Y":
                        rel_final = "Y"
            if (rel_final == 'N'):
                df.loc[(df.seq_id == row.seq_id) & (~(row.seq_id in fps)) & (~(row.seq_id in fns)), 'TN'] = 1
        tns = np.unique(df.seq_id[(df.TN == 1)])
        tns = list(set(tns))
        num_tn = len(tns)
    
    num_fusions = len(np.unique(df.seq_id[(df['genomic_subtype_std_rel'].str.contains("Y")) | (df['gu_fusion_std_rel'].str.contains("Y")) | (df['gu_pax5_status_std_rel'].str.contains("Y"))]))
    
    # check length of total matches metadata total
    if metadata_fusion == True:
        total = tps + fns + tns + fps
        total = set(total)
        assert(len(total) == u1), "ERROR: number of relevant seq_ids in metadata and evaluated seq_ids do not match!"
        assert(set(tps).isdisjoint(set(fns))), "ERROR: overlap between TPs and FNs"
        assert(set(tps).isdisjoint(set(tns))), "ERROR: overlap between TPs and TNs"
        assert(set(tps).isdisjoint(set(fps))), "ERROR: overlap between TPs and FPs"
        assert(set(fns).isdisjoint(set(tns))), "ERROR: overlap between FNs and TNs"
        assert(set(fns).isdisjoint(set(fps))), "ERROR: overlap between FNs and FPs"
        assert(set(tns).isdisjoint(set(fps))), "ERROR: overlap between TNs and FPs"
    

        print("---------------------------")
        print(f"{name} performance")
        print("---------------------------")
        print(f"Initial number of seq_ids: {num_seqids}")
        print(f"Number of seq_ids filtered because lineage was not in {ls}: {l0-l1}")
        print(f"Number of seq_ids with unknown genomic_std_subtypes: {u0-u1}")
        print(f"Number of seq_ids after filtering: {u1}")
        print(f"Number of seq_ids with relevant fusions (should match TPs + FNs): {num_fusions}")
        print(f"True Positives: {num_tp}")
        print(f"False Negatives: {num_fn}")
        print(f"False Positives: {num_fp}")
        print(f"True Negative: {num_tn}")
        print(f"Sanity Check length of total (should be num of seq_ids after filtering): {len(total)}")
        if num_fn+num_tp == 0:
            se = np.nan
        else:
            se = num_tp/(num_fn+num_tp)
        if num_tn+num_fp == 0:
            sp = np.nan
        else:
            sp = num_tn/(num_tn+num_fp)
        if num_tp+num_fp == 0:
            pr = np.nan
        else:
            pr = num_tp/(num_tp+num_fp)
        if 2*num_tp+num_fp+num_fn == 0:
            f1 = np.nan
        else:
            f1 = 2*num_tp/(2*num_tp+num_fp+num_fn)
        print(f"Sensitivity (Recall): {se:.4f}")
        print(f"Specificity: {sp:.4f}")
        print(f"Precision: {pr:.4f}")
        print(f"F1 stat: {f1:.4f}")
    
    
    df.to_csv(fd_output_path+'/'+name+'.csv', sep="\t")
    
    if metadata_fusion == True:
        return tps, fns, tns, fps, se, sp, pr, f1

# function to plot subplots of fusion results
def plot_results(metadata_fusion, fd_output_path, name, label_size=1, xtick_size=0.05, xtick_gap=0.5, xtick_rot=45, ytick_size=1.5, ytick_gap=0.5):
    df = pd.read_csv(fd_output_path+'/'+name+'.csv', sep="\t")
    if metadata_fusion == True:
        tps = list(set(np.unique(df.seq_id[(df['ges_fusion_yn'] == "Y") & (df.TP_seq_id == 1)])))
        fns = list(set(np.unique(df.seq_id[df.FN == 1])))
        fps = list(set(np.unique(df.seq_id[df['FP_seq_id'] == 1])))
        # which one of these from the set above are NOT in false positives...these are truly negative on seq_id level
        tns = np.unique(df.seq_id[df['TN'] == 1])
        tns = list(set(tns).difference(fps))
        # plot the barchart num reads distribution of fusions detected for each sample
        plt.figure(figsize=(1000, 600))
        rel_seq_ids = sorted(np.append(tps, fns))
        rel_seq_ids = sorted(np.append(rel_seq_ids, fps))
    else:
        rel_seq_ids = sorted(list(set(np.unique(df.seq_id))))
    # try to create nxn matrix of plots
    numcols = math.ceil(math.sqrt(len(rel_seq_ids)))
    fig, ax = plt.subplots(nrows=numcols, ncols=numcols)
    row_index = 0
    col_index = 0
    plt.subplots_adjust(hspace=1)
    fig.suptitle(name)
    for seq_id in rel_seq_ids:
        temp = df.loc[(df['seq_id'] == seq_id),:].copy()
        # highlight the true positives
        colors = ['blue'] * temp.shape[0]
        if metadata_fusion == True:
            js = np.where(temp.TP_fusion == 1)[0]
            ks = np.where((temp.FP_fusion == 1) & (temp.fusion_detected_rel == 'Y'))[0]
            if len(js) > 0:
                for j in js:
                    colors[j] = 'red'
            if len(ks) > 0:
                for k in ks:
                    colors[k] = 'purple'
        ax[row_index][col_index].bar(x='fusion_detected_std12', height='num_reads', color=colors, data=temp)
        for xtick, color in zip(ax[row_index][col_index].get_xticklabels(), colors):
            xtick.set_color(color)
        ax[row_index][col_index].xaxis.set_tick_params(labelsize=xtick_size, labelrotation=xtick_rot)
        ax[row_index][col_index].yaxis.set_tick_params(labelsize=ytick_size, length=1)
        ax[row_index][col_index].tick_params(bottom=False) 
        ax[row_index][col_index].tick_params(axis="x",direction="in", pad=xtick_gap)
        ax[row_index][col_index].tick_params(axis="y",direction="in", pad=ytick_gap)
        ax[row_index][col_index].text(0.125, 0.5, seq_id, fontsize=label_size, transform=ax[row_index][col_index].transAxes)
        # plt.ylabel('num_reads')
        plt.show()
        # if df[(df['seq_id'] == seq_id) & (df.TP == 1)].shape[0] > 0:
            # df[(df['seq_id'] == seq_id) & (df.TP == 1)].plot.bar(x='fusion_detected_std12', y='num_reads', color='red')
        # plt.show()
        if col_index == numcols-1:
            col_index = 0
            row_index += 1
        else:
            col_index += 1
    fig.savefig(fd_output_path+"/"+name+".png", dpi=192*5, facecolor="white")  
    
# function to determine strand given alignment:
def det_strand(aln):
    if aln.rv == True:
        strand = '-'
    else:
        strand = '+'
    return strand

    
# function to determine breakpoints...given two alignments (needs paf module)
# returns the breakpoitn of whatever is listed first ie aln0
def det_bp(aln0, aln1):
    # assert(aln0.qs != aln1.qs), "Alignments have equal start positions. We cannot discriminate which is first..."
    # first determine which alignment is first on the read ie reorder if necessary
    if aln0.qs == aln1.qs:
        return np.nan
    elif aln0.qs < aln1.qs:
        temp0 = aln0
        temp1 = aln1
    else:
        temp0 = aln1
        temp1 = aln0
        # temp0 = aln0
        # aln0 = aln1
        # aln1 = temp
    # define the strands
    temp0_s = det_strand(temp0)
    temp1_s = det_strand(temp1)
    # get the strand pair to later define the breakpoints
    orts = temp0_s+temp1_s
    if orts == '++':
        bp0 = temp0.te
        bp1 = temp1.ts
    elif orts == '+-':
        bp0 = temp0.te
        bp1 = temp1.te
    elif orts == '-+':
        bp0 = temp0.ts
        bp1 = temp1.ts
    elif orts == '--':
        bp0 = temp0.ts
        bp1 = temp1.te
    if aln0 == temp0:
        return bp0
    elif aln0 == temp1:
        return bp1

class CandGene:
    def __init__(self, name, chrom, start, end):
        self.name = name
        self.chr = chrom
        self.s = start
        self.e = end
    def __hash__(self):
        return hash(self.name)

# distance between alignments on query
def aln_dist(h0, h1):
    return min(h0.qe, h1.qe) - max(h0.qs, h1.qs)

# distance between targets if on same chromosome
def t_dist(h0, h1):
    assert(h0.t == h1.t), "Alignments need to be on same chromosome for this metric!"
    return min(h0.te, h1.te) - max(h0.ts, h1.ts)

# percentage of overlapping alignments. as a percentage of aln1 so order of args matters. returns percentage as xx.xx%
def q_overlap(aln0, aln1):
    return (min(aln0.qe, aln1.qe) - max(aln0.qs, aln1.qs))/(aln1.qe - aln1.qs) * 100
    

# overlap between gene and alignment
def gh_overlap(cgene_s, cgene_e, aln): # ie gene start, gene end, alignment paf record
    return min(cgene_e, aln.te) - max(cgene_s, aln.ts)

# create a dataframe of supporting reads and read stats for your cohort
# rdsupport_path is the path to directory containing pkl objects of supporting reads for fusion caller output
def agg_support(metadata_path, fusion_std_path, gu_std_fusion_path, gu_pax5_status_std_path, 
                rdsupport_path,
                min_anch=50, max_gap=30, max_overlap=50, min_gene_gap=10, 
                max_gene_overlap=10, qmax_overlap=40, bp_win=10, min_ct=2):
    df = read_md(metadata_path, fusion_std_path, gu_std_fusion_path, gu_pax5_status_std_path)
    # set the relevant params for determining filter success
    # bin_size=1000000
    min_anchor=50
    max_gap=30
    max_overlap=50
    min_gene_gap=10
    max_gene_overlap=10
    qmax_overlap = 40
    bp_win=10
    min_ct=2
    # wsize = 50
    # hwsize = wsize/2

    # initialize lists
    seq_met = defaultdict(list)
    seq_id_df = []
    qid_df = []
    g0_df = []
    g1_df = []
    h0l_df = []
    h1l_df = []
    a0_df = []
    a1_df = []
    h0h1d_df = []
    t0t1d_df = []
    q0ovlp_df = []
    q1ovlp_df = []
    bp0d_df = []
    bp1d_df = []
    a0_pf = []
    a1_pf = []
    h0h1d_pf = []
    t0t1d_pf = []
    q0ovlp_pf = []
    q1ovlp_pf = []
    bp0d_pf = []
    bp1d_pf = []

    for seq_id in list(np.unique(df.seq_id)):
        # load the objects
        with open(rdsupport_path+"/"+seq_id+'_data.pkl', 'rb') as inp:
            fg_lookup = pickle.load(inp)
            bp_lookup = pickle.load(inp)
            fr_lookup = pickle.load(inp)
            win_lookup = pickle.load(inp)
            fusions = pickle.load(inp)

        for f in fg_lookup:
            seq_met[f] = []
            # genes
            cg0 = fg_lookup[f][0]
            cg1 = fg_lookup[f][1]
            g0 = cg0.name
            g1 = cg1.name
            g_list = [g0, g1]
            cg_list = [cg0, cg1]
            # alignments
            a_list = []
            m_list = []
            for algn_pair in fr_lookup[f]:
                seq_id_df.append(seq_id)
                h0 = algn_pair[0]
                h1 = algn_pair[1]
                qid_df.append(h0.q)
                h0l_df.append(h0.ql)
                h1l_df.append(h1.ql)
                a_list.append([h0, h1])
                # metrics
                # overlap between gene and alignment
                a0 = gh_overlap(cg0.s, cg0.e, h0)
                a1 = gh_overlap(cg1.s, cg1.e, h1)
                # distance between alignments on query
                h0h1d = aln_dist(h0, h1)
                # distance between targets if on same chr
                if h0.t == h1.t:
                    t0t1d = t_dist(h0, h1)
                else:
                    t0t1d =  np.nan
                # percentage of overlapping alignments with respect to given alignment length
                q0ovlp = q_overlap(h1, h0)
                q1ovlp = q_overlap(h0, h1)
                # distance between breakpoint and nearest gene boundary
                bp0 = det_bp(h0, h1)
                bp1 = det_bp(h1, h0)
                # if inside the gene return np.nan...since not an issue. 
                # if outside the gene...return the distance
                if (bp0 >= cg0.s) and  (bp0 <= cg0.e):
                    bp0d = np.nan
                else:
                    if bp0 < cg0.s:
                        bp0d = cg0.s - bp0
                    elif bp0 > cg0.e:
                        bp0d = bp0 - cg0.e
                if (bp1 >= cg1.s) and  (bp1 <= cg1.e):
                    bp1d = np.nan
                else:
                    if bp1 < cg1.s:
                        bp1d = cg1.s - bp1
                    elif bp1 > cg1.e:
                        bp1d = bp1 - cg1.e
                m_list.append([a0, a1, h0h1d, t0t1d, q0ovlp, q1ovlp, bp0d, bp1d])

                g0_df.append(g0)
                g1_df.append(g1)

                a0_df.append(a0)
                if a0 > min_anchor:
                    a0_pf.append(True)
                else:
                    a0_pf.append(False)

                a1_df.append(a1)
                if a1 > min_anchor:
                    a1_pf.append(True)
                else:
                    a1_pf.append(False)

                h0h1d_df.append(h0h1d)
                if h0h1d > max_overlap or h0h1d < -max_gap:
                    h0h1d_pf.append(False)
                else:
                    h0h1d_pf.append(True)

                t0t1d_df.append(t0t1d)
                if t0t1d > -min_gene_gap or t0t1d > max_gene_overlap:
                    t0t1d_pf.append(False)
                else:
                    t0t1d_pf.append(True)

                q0ovlp_df.append(q0ovlp)
                if q0ovlp > qmax_overlap:
                    q0ovlp_pf.append(False)
                else:
                    q0ovlp_pf.append(True)
                q1ovlp_df.append(q1ovlp)
                if q1ovlp > qmax_overlap:
                    q1ovlp_pf.append(False)
                else:
                    q1ovlp_pf.append(True)
                bp0d_df.append(bp0d)
                bp1d_df.append(bp1d)
                if bp0 < (cg0.s - bp_win) or bp0 > (cg0.e + bp_win):
                    bp0d_pf.append(False)
                else:
                    bp0d_pf.append(True)
                if bp1 < (cg1.s - bp_win) or bp1 > (cg1.e + bp_win):
                    bp1d_pf.append(False)
                else:
                    bp1d_pf.append(True)


        seq_met[f] = [g_list, a_list, m_list]
        results = pd.DataFrame({
              'seq_id': seq_id_df,
              'qid': qid_df,
              'g0': g0_df, 
              'g1': g1_df, 
              'h0_qlen': h0l_df, 
              'h1_qlen': h1l_df, 
              'g0_h0_ovlp': a0_df, 
              'g0_h0_ovlp_pf': a0_pf,
              'g1_h1_ovlp': a1_df,
              'g1_h1_ovlp_pf': a1_pf,
              'h0h1_dist': h0h1d_df,
              'h0h1_dist_pf': h0h1d_pf,
              't0t1_dist': t0t1d_df,
              't0t1_dist_pf': t0t1d_pf,
              'q0_ovlp_perc': q0ovlp_df,
              'q0_ovlp_perc_pf': q0ovlp_pf,
              'q1_ovlp_perc': q1ovlp_df,
              'q1_ovlp_perc_pf': q1ovlp_pf,
              'bp0_dist': bp0d_df,
              'bp0_dist_pf': bp0d_pf,
              'bp1_dist': bp1d_df,
              'bp1_dist_pf': bp1d_pf
             })
        
    results.loc[:, 'overall_pf'] = results.g0_h0_ovlp_pf & results.g1_h1_ovlp_pf & results.h0h1_dist_pf & results.t0t1_dist_pf & results.q0_ovlp_perc_pf & results.q1_ovlp_perc_pf & results.bp0_dist_pf & results.bp1_dist_pf
    return results



# investigate supporting reads for a fusion
def read_lkup(results, seq_id, gene1, gene2):
    o = results.loc[(results.seq_id == seq_id) & 
                (((results.g0 == gene1) & (results.g1 == gene2))|
                 ((results.g1 == gene1) & (results.g0 == gene2))), :]
    return o

