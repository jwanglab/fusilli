import sys
sys.path.append('/nas/longleaf/home/jclin/jwanglab_users_jclin/nanopore-dx/10_fusion_detection/fusion_detection_modules')
import common as c
import pandas as pd
import numpy as np

# inputs
fd_output_path = sys.argv[1]
fusion_master_path = sys.argv[2]
fusion_std_path = sys.argv[3]
name = sys.argv[4]
metadata_path = sys.argv[5]
metadata_fusion = sys.argv[6]
bin_size = int(sys.argv[7])
min_anchor = int(sys.argv[8])
max_gap = int(sys.argv[9])
max_overlap = int(sys.argv[10])
min_gene_gap = int(sys.argv[11])
max_gene_overlap = int(sys.argv[12])
qmax_overlap = int(sys.argv[13])/100
bp_win = int(sys.argv[14])
min_ct = int(sys.argv[15])
cohort = sys.argv[16]
filt = sys.argv[17]
comment = sys.argv[18]
gu_fusion_std_path = sys.argv[19]
gu_pax5_status_std_path = sys.argv[20]
wlist_path = sys.argv[21]


if metadata_fusion == "True":
    metadata_fusion = True
elif metadata_fusion == "False":
    metadata_fusion = False
if filt == "True":
    filt = True
elif filt == "False":
    filt = False


    
    
# read in metadata and do some standardization of fusion names
metadata = c.read_md(metadata_path, fusion_std_path, gu_fusion_std_path, gu_pax5_status_std_path)

# aggregate results from as_karyo_fusions into a dataframe
seq_ids = np.asarray(metadata.seq_id)
gene1_list = []
gene2_list = []
num_reads_list = []
seq_id_list = []
for seq_id in seq_ids:
    i = 0
    start = False
    with open(fd_output_path + '/' + seq_id + '_fusilli.txt') as fp:
        for line in fp:
            if line.strip() == 'No fusions detected!':
                print(f"{seq_id} has no fusions called")
                gene1_list.append(np.nan)
                gene2_list.append(np.nan)
                num_reads_list.append(np.nan)
                seq_id_list.append(seq_id)
                break # go to next seq_id
            else:
                gene1 = line.strip().split("\t")[0]
                gene2 = line.strip().split("\t")[1]
                num_reads = line.strip().split("\t")[2]
                gene1_list.append(gene1)
                gene2_list.append(gene2)
                num_reads_list.append(num_reads)
                seq_id_list.append(seq_id)

df = pd.DataFrame({'seq_id': seq_id_list,
             'gene1': gene1_list,
             'gene2': gene2_list,
             'num_reads' : num_reads_list}
            )

if df.shape[0] == 0:
    print("There were no fusions identified in this cohort. Stopping aggregation.")
    sys.exit()

df['fusion_detected_std12'] =  df.gene1.apply(str) + "::" + df.gene2.apply(str)
df['fusion_detected_std21'] =  df.gene2.apply(str) + "::" + df.gene1.apply(str)

# read in gene fusion master
fm = c.read_fm(fusion_master_path)
# compare fusions in fusion_std against fusion_master
c.compare_fstd_fm(fusion_std_path, fm)
# compare fusions in metadata against fusion_master
c.compare_md_fm(metadata, fm)


# determine overall performance
tps, fns, tns, fps, se, sp, pr, f1 = c.def_results(df, metadata, fm, name, fd_output_path, metadata_fusion=metadata_fusion, fusion_std_path=fusion_std_path, gu_fusion_std_path=gu_fusion_std_path, gu_pax5_status_std_path=gu_pax5_status_std_path, wlist_path=wlist_path)
# plot the results
c.plot_results(metadata_fusion, fd_output_path, name)

# write the results
with open('/nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/07_new_fd/result_summary.tsv', 'a') as f:
    f.write(cohort + "\t" + str(se) + "\t" + str(sp) + "\t" + str(pr) + "\t" + str(f1) + "\t" + str(bin_size) + "\t" + str(min_anchor) + "\t" + str(max_gap) + "\t" + str(max_overlap) + "\t" + str(min_gene_gap) + "\t" + str(max_gene_overlap) + "\t" + str(qmax_overlap) + "\t" + str(bp_win) + "\t" + str(min_ct) + "\t" + str(filt) + "\t" + str(comment) + "\n")