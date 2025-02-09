from collections import defaultdict
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import collections
import pickle
import paf
import ga_calc
import gene
import sys
from tqdm import tqdm
import time
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

chroms = ["NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12", "NC_000005.10", "NC_000006.12", 
          "NC_000007.14", "NC_000008.11", "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12",
          "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10", "NC_000017.11", "NC_000018.10", 
          "NC_000019.10", "NC_000020.11", "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10"]
chr_idx = {chroms[c]:c for c in range(len(chroms))}

chr_len = [None for c in chroms]
chr_bins = [None for c in chroms]

# reads bed file and outputs dictionary of genes with genomic positions
def read_bed(bed_file):
    genes = {}
    for line in open(bed_file):
        chrom, st, en, gene = line.strip().split('\t')
        if chrom not in chroms:
            continue
        genes[gene] = (chroms.index(chrom), int(st), int(en))
    return genes

def read_fmaster(fusion_master_file):
    result = []
    for line in open(fusion_master_file):
        gene1, gene2 = line.strip().split('\t')
        result.append(gene1 + "::" + gene2)
    return result


def main(paf_file, bed_file, fm_file, nfusm, outpath, min_anchor, max_gap, max_overlap, max_gene_overlap, bp_win, qmax_overlap, filt, min_ct, rep):
    
    genes = read_bed(bed_file)
    fm_list = read_fmaster(fm_file)
    hits = defaultdict(list)
    sample_id=paf_file.split('/')[-1].split('.paf')[0]
    if outpath != '':
        outpath1 = outpath + '/' + sample_id + '_fusilli.txt'
    with tqdm(total=sum(1 for _ in open(paf_file)), desc="analyzing the PAF file at " + paf_file + "...") as pbar:
        with open(paf_file) as pf:
            for line in pf:
                # create aln object from each line in the paf file
                try:
                    a = paf.aln(line)
                except:
                    print("WARNING: PAF file has bad format or is truncated, stopped reading it here!")
                    break
                # if the target sequence name is in the list of chromosomes we are interested in...
                if a.t in chr_idx:
                    c = chr_idx[a.t]
                    for g in genes:
                        #...and the gene overlaps with the alignment by at least min_anchor...
                        if genes[g][0] == c and ga_calc.ga_ovlp(genes[g][1], genes[g][2], a) > min_anchor:
                            #...add the alignment to the hits dictionary with the query sequence hash as the key
                            if a.q in hits:
                                hits[a.q].append(a)
                            else:
                                # hits is a default dict of form hits = [(a.q, [alignments])] where a.q is the key ie query sequence hash followed by list of alignments
                                hits[a.q] = [a]
                pbar.update(1)
    # fusions is a default dict of form fusions = [(a, b)] where a and b are the gene names of the fusion
    # fusions = defaultdict(int)
    fr_lookup = defaultdict(list)

    rs = []
    rls = []
    gene0l = []
    gene1l = []
    chr0s = []
    chr1s = []
    f1 = []
    f2 = []
    f3 = []
    f4 = []
    f5 = []
    f6 = []
    f1b = []
    f2b = []
    f3b = []
    f4b = []
    f5b = []
    f6b = []
    f7b = []
    fus = []
    # fo = []
    for read in hits: 
        # a single query ie read can have multiple alignments or mappings in the genome
        # hits are reads where they fall within a target chromosome
        for i in range(len(hits[read])): # len of hits[read] is the number of alignments associated with given query seq hash/read
            a0 = hits[read][i] # a0 is the ith or first alignment...this is basically a line from the PAF file
            c0 = chr_idx[a0.t] # a0.t is the corresponding target sequence name ie chromosome, so c0 is the index of the chr...chr1 --> 0, chr2 --> 1 etc
            for j in range(i+1, len(hits[read])): # for the NEXT (ie i+1) read associated with the given query....also get the alignment and chr index
                a1 = hits[read][j]
                c1 = chr_idx[a1.t]
                # get the length of the overlap between the two alignments....are they greater than 100? if yes, then they overlap
                # qe is query end and qs is query start with respect to the read! so a read starts at 0 ... length(read)
                
                '''
                # filter - appropriate gap in the alignment records
                if filt:
                    if ga_calc.a_dist_tf(a0, a1, max_overlap, max_gap) == False:
                        continue
                # filter - genes must not overlap too much
                if filt:
                    if c0 == c1 and (ga_calc.t_dist(a0, a1) > max_gene_overlap):
                        continue
                # filter - overlapping portion must not be more than 40% of query lengths
                if filt:
                    if (ga_calc.q_overlap(a1, a0) > qmax_overlap) or (ga_calc.q_overlap(a0, a1) > qmax_overlap):
                        continue
                '''
                
                # these are generic identifiers of the genes involved based on chr and halfway point location
                # g0 = (c0, (a0.ts+(a0.te-a0.ts)//2)//1000000) # (chrom, loc[Mbp]) this is a tuple...representing the chrom and halfway point of each alignment??
                # g1 = (c1, (a1.ts+(a1.te-a1.ts)//2)//1000000) 
                g0 = None
                g1 = None
                # the actual gene names are assigned further down below

                # assign gene names instead of coordinates to those matching one of our known targets
                
                # check if there is more than one gene for each alignment...??
                g0l = []
                g1l = []
                
                for g in genes:
                    # if the chr is the same...and the overlap between the gene and target sequence is greater than min_anchor...assign g0 to given g ie of class string
                    

                        
                    if genes[g][0] == c0 and ga_calc.ga_ovlp(genes[g][1], genes[g][2], a0) > min_anchor: ## ?? this may not necessarily cover the breakpoint...ex if aln is entirely within gene
                        g0l.append(g)
                        # g0 = g
                        g0_c = genes[g][0]
                        g0_s = genes[g][1]
                        g0_e = genes[g][2]
                        cg0 = gene.CandGene(g0, g0_c, g0_s, g0_e)
                    if genes[g][0] == c1 and ga_calc.ga_ovlp(genes[g][1], genes[g][2], a1) > min_anchor:
                        g1l.append(g)
                        # g1 = g
                        g1_c = genes[g][0]
                        g1_s = genes[g][1]
                        g1_e = genes[g][2]
                        cg1 = gene.CandGene(g1, g1_c, g1_s, g1_e)
                # print(g0, g1)
                # if len(set(g0l)) > 1:
                #     g0 = g0l
                #     # print(g0l)
                # if len(set(g1l)) > 1:
                #     g1 = g1l
                #     print(g1l)
                
                
                '''
                # breakpoint must be within 10 bp of gene start/end
                if filt:
                    if (ga_calc.det_bp(a0, a1) < (g0_s - bp_win)) | (ga_calc.det_bp(a0, a1) > (g0_e + bp_win)):
                        continue
                    if (ga_calc.det_bp(a1, a0) < (g1_s - bp_win)) | (ga_calc.det_bp(a1, a0) > (g1_e + bp_win)):
                        continue
                '''
                # only consider genes on different alignments...an alignment can have multiple genes...but here we only consider gene pairs from different alignments ie g0 from a0 and g1 from a1
                # ?? this is a limitation of the software...the only way to get metrics for multiple genes from the same target alignment is to rerun the mapping with a reference broken down by gene target starts and ends
                for k in range(len(g0l)): 
                    g0 = g0l[k]
                    for l in range(len(g1l)):
                        g1 = g1l[l]
                        # only consider genes when they are different partner genes
                        if g0 == g1 or ((g0 == None) or (g1 == None)):
                            continue
                        key = (g0, g1)
                        r_val = (a0, a1)
                        # standardize the fusion name here
                        if (g0 < g1):
                        #     # fus0.append(g0)
                        #     # fus1.append(g1)
                            key = (g0, g1)

                        #     # g_val = [cg0, cg1]
                            r_val = (a0, a1)
                        else:
                        #     # fus0.append(g1)
                        #     # fus1.append(g0)
                            key = (g1, g0)
                        #     # g_val = [cg1, cg0]
                            r_val = (a1, a0)
                        # if the read has already been used in a supporting alignment pair, skip, since we are counting only the number of supporting reads
                        # if (len(fr_lookup[key])!= 0):
                        #     if r_val[0].q == fr_lookup[key][0][0].q:
                        #         continue
                        # fg_lookup[key] =  g_val
                        fr_lookup[key].append(r_val)
                        # win_lookup[key] = defaultdict(list)
                        # fusions[key] += 1
                        fus.append(key)
                        gene0l.append(g0)
                        gene1l.append(g1)
                        chr0s.append(g0_c)
                        chr1s.append(g1_c)
                        rs.append(a0.q)
                        rls.append(a0.ql)
                        f1.append(ga_calc.a_dist(a0, a1))
                        if a0.t == a1.t:
                            f2.append(ga_calc.t_dist(a0, a1))
                        else:
                            f2.append(None)
                        f3.append(ga_calc.q_overlap(a0, a1))
                        f4.append(ga_calc.q_overlap(a1, a0))
                        f5.append(ga_calc.det_bp(a0, a1))
                        f6.append(ga_calc.det_bp(a1, a0))

                        b1 = ga_calc.a_dist_tf(a0, a1, max_overlap, max_gap)
                        b2 = ga_calc.t_dist_tf(a0, a1, max_gene_overlap)
                        b3 = ga_calc.q_overlap_tf(a0, a1, qmax_overlap)
                        b4 = ga_calc.q_overlap_tf(a1, a0, qmax_overlap)
                        b5 = ga_calc.bp_tf(a0, a1, g0_s, g0_e, bp_win)
                        b6 = ga_calc.bp_tf(a1, a0, g1_s, g1_e, bp_win)

                        if nfusm:
                            if ((key[0] + "::" + key[1] in fm_list) or (key[1] + "::" + key[0] in fm_list)):
                                b7 = True
                            else:
                                b7 = False
                        else:
                            b7 = True
                        
                        f1b.append(b1)
                        f2b.append(b2)
                        f3b.append(b3)
                        f4b.append(b4)
                        f5b.append(b5)
                        f6b.append(b6)
                        f7b.append(b7)
                        # f7b.append(b7)
                        # fo.append(bo)


    # fusions = sorted([(f,fusions[f]) for f in fusions], key=lambda a:a[1])
    # fusions_filt = {}

    results = pd.DataFrame({
              'gene0': gene0l,
              'gene1': gene1l,
              'chr0': chr0s,
              'chr1': chr1s,
              'fusion': fus,
              'read_query': rs,
              'read_length': rls,
              'aln_dist_rd': f1, 
              'aln_dist_rd_tf': f1b, 
              'target_dist': f2,
              'target_dist_tf': f2b, 
              'ovlp_perc_aln0': f3, 
              'ovlp_perc_aln0_tf': f3b, 
              'ovlp_perc_aln1': f4, 
              'ovlp_perc_aln1_tf': f4b, 
              'aln0_bp': f5, 
              'aln0_bp_win_tf': f5b, 
              'aln1_bp': f6, 
              'aln1_bp_win_tf': f6b,
              'in_fus_mast_tf': f7b,
             })
    
    # boolean for overall filter pass True/False
    if filt:
        results['overall_filt_tf'] = np.where((results['aln_dist_rd_tf'] == True) & (results['target_dist_tf'] == True) & (results['ovlp_perc_aln0_tf'] == True) & (results['ovlp_perc_aln1_tf'] == True) & (results['aln0_bp_win_tf'] == True) & (results['aln1_bp_win_tf'] == True), True, False)
    else:
        results['overall_filt_tf'] = True
    
    # count the number of distinct passing-all-filters reads supporting each fusion 
    ct_lkup = results.loc[results.overall_filt_tf == True, ['fusion', 'read_query']].groupby('fusion')['read_query'].nunique().reset_index(name='ct')
    results['ct'] = results['fusion'].map(ct_lkup.set_index('fusion')['ct'])
    results['ct_tf'] = results['ct'].apply(lambda x: x >= min_ct)
    results['overall_tf'] = np.where((results['overall_filt_tf'] == True) & (results['in_fus_mast_tf'] == True) & (results['ct_tf'] == True), True, False)
    # results['overall_tf'] = results['overall_filt_tf'] == results['in_fus_mast_tf'] == results['ct_tf'] == True


    if outpath != '':
        orig_stdout = sys.stdout
        fo = open(outpath1, 'w')
        sys.stdout = fo
        if results.overall_tf.any():
            results_sub = results.loc[results.overall_tf == True, ['fusion', 'ct']].drop_duplicates().sort_values(by='ct', ascending=False)
            for index, row in results_sub.iterrows():
                print(row.fusion[0] + "\t" + row.fusion[1] + "\t" + str(row.ct))
        else:
            print("No fusions detected!")
        if rep:
            results.to_csv(outpath + '/' + sample_id + '_fusilli_read_summary.txt', sep="\t", index=False)
        sys.stdout = orig_stdout
        fo.close()

    print()
    print("---------------------------")
    print(" FUSILLI Fusion Detection")
    print("---------------------------")
    if results.overall_tf.any():
        results_sub = results.loc[results.overall_tf == True, ['fusion', 'ct']].drop_duplicates().sort_values(by='ct', ascending=False)
        for index, row in results_sub.iterrows():
            print(row.fusion, row.ct)
    else:
        print("No fusions detected!")
    if rep:
        print("---------------------------")
        print(" FUSILLI Fusion Read Summary")
        print("---------------------------")
        print(results)
    print()
    print("finished!")

    # print(results)

    # for k, ct in fusions:
    #     if ct >= min_ct:
    #         if nfusm:
    #             if ((k[0] + "::" + k[1] in fm_list) or (k[1] + "::" + k[0] in fm_list)):
    #                 fusions_filt[k] = ct
    #         else:
    #             fusions_filt[k] = ct
    # for index, row in results.iterrows():
    #     if row.fo == True:
    #         print(row)
    


    


    # if outpath != '':
    #     orig_stdout = sys.stdout
    #     fo = open(outpath, 'w')
    #     sys.stdout = fo
    #     print()
    #     print("---------------------------")
    #     print(" Fusion detection")
    #     print("---------------------------")
    #     if len(fusions_filt) > 0:
    #         for f,ct in fusions_filt.items():
    #             print(str(f), ct)
    #     else:
    #         print("No fusions detected!")
    #     sys.stdout = orig_stdout
    #     fo.close()
    # else:
    #     print()
    #     print("---------------------------")
    #     print(" Fusion detection")
    #     print("---------------------------")
    #     if len(fusions_filt) > 0:
    #         for f,ct in fusions_filt.items():
    #             print(f, ct)
    #     else:
    #         print("No fusions detected!")
    
    # print()
    # print("finished!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser("B-ALL fusion caller based on ONT RNA-seq", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--paf", help="PAF of reads to hg38 (mm2 -x map-ont)")

    parser.add_argument("-b", "--bed", help="BED file including target genes, default is BED file specific for B-ALL", default='/proj/jwanglab/users/jclin/nanopore-dx/10_fusion_detection/fusilli/data/peds_leuk_genes.240422.bed')
    parser.add_argument("-fm", "--fusionmaster", help="tab-delimited file with fusion genes of interest, default is file for B-ALL fusions", default='/proj/jwanglab/users/jclin/nanopore-dx/10_fusion_detection/fusilli/data/ball_fusion_master.txt')
    parser.add_argument("-nfm", "--no_fusion_master", action='store_false', help="if included argument, any gene pairs from the bed file are considered")
    parser.add_argument("-o", "--outpath", help="if included argument, output path for fusion detection output", default='')
    parser.add_argument("-a", "--anchor", type=int, help="minimum anchor size ie overlap between read and target gene", default=50)
    parser.add_argument("-mxg", "--maxgap", type=int, help="maximum gap (bp) between alignments", default=30)
    parser.add_argument("-mxo", "--maxoverlap", type=int, help="maximum overlap (bp) between alignments", default=50)
    parser.add_argument("-mxgo", "--maxgeneoverlap", type=int, help="maximum overlap (bp) between genes if target genes on same chromosome", default=10)
    parser.add_argument("-bpw", "--bpwin", type=int, help="breakpoint must be within this number of bases before start and end of gene", default=10)
    parser.add_argument("-qmo", "--qmaxoverlap", type=int, help="overlapping portion must not be more than certain percentage of a given read length", default=40)
    parser.add_argument("-nf", "--no_filt", action='store_false', help="if included argument, no filters applied, (except for mincount)")
    parser.add_argument("-mnc", "--mincount", type=int, help="minimum number reads to count as a fusion", default=2)
    parser.add_argument("-r", "--report",  action='store_true', help="if included argument, a report will be generated of all reads, what filters passed and whether they are fusions not included in a fusion master")
    
    args = parser.parse_args()

    main(args.paf, 
         args.bed,
         args.fusionmaster,
         args.no_fusion_master,
         args.outpath,
         args.anchor, 
         args.maxgap, 
         args.maxoverlap, 
         args.maxgeneoverlap, 
         args.bpwin,
         args.qmaxoverlap,
         args.no_filt,
         args.mincount,
         args.report
         )