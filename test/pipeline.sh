#!/bin/bash

module load anaconda
conda activate fusion

input_path=/nas/longleaf/home/jclin/jwanglab_users_jclin/data/xx_ncbi_grch38_genome_alignments
# old /proj/jwanglab/projects/nanopore_dx/dna_adaptive/peds_leuk_genes.bed
# from brady and klco /proj/jwanglab/projects/nanopore_dx/dna_adaptive/peds_leuk_genes.240422.bed
bed_path=/proj/jwanglab/projects/nanopore_dx/dna_adaptive/peds_leuk_genes.240422.bed
# dux4 mods /proj/jwanglab/projects/nanopore_dx/dna_adaptive/peds_leuk_genes.240515.bed
# bed_path=/nas/longleaf/home/jclin/jwanglab_users_jclin/data/references/ncbi_grch38/GCF_000001405.40/genomic_subset.bed
# bed_path=/proj/jwanglab/projects/nanopore_dx/dna_adaptive/peds_leuk_genes.240827.txt
fusion_master_path=/nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/fusion_master_20240415.txt
fusion_std_path=/nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/fusion_std_20240613.txt
gu_fusion_std_path=/nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/gu_fusion_std_20240620.txt
gu_pax5_status_std_path=/nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/gu_pax5_status_std_20240620.txt
wlist_path1=none
wlist_path2=/nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/whitelist.txt
# low depth is metadata_low_depth.txt
# high depth is metadata_high_depth.txt
metadata_path=/proj/jwanglab/users/jclin/data/metadata_high_depth.txt
metadata_fusion=True
seq_id_col=0
cohort=high_depth_fusilli
comment='250208_fusilli_test_high_depth'

bin_size_arr=(1000000)
min_anchor_arr=(50)
max_gap_arr=(30)
max_overlap_arr=(50)
min_gene_gap_arr=(10)
max_gene_overlap_arr=(10)
qmax_overlap_arr=(40)
bp_win_arr=(10)
min_ct_arr=(2)
filt=True

len=${#bin_size_arr[@]}; echo $len
len=${#min_anchor_arr[@]}; echo $len
len=${#max_gap_arr[@]}; echo $len
len=${#max_overlap_arr[@]}; echo $len
len=${#min_gene_gap_arr[@]}; echo $len
len=${#max_gene_overlap_arr[@]}; echo $len
len=${#qmax_overlap_arr[@]}; echo $len
len=${#bp_win_arr[@]}; echo $len
len=${#min_ct_arr[@]}; echo $len



for (( i=0; i<${len}; i++ ));
do
    bin_size=${bin_size_arr[$i]}
    min_anchor=${min_anchor_arr[$i]}
    max_gap=${max_gap_arr[$i]}
    max_overlap=${max_overlap_arr[$i]}
    min_gene_gap=${min_gene_gap_arr[$i]}
    max_gene_overlap=${max_gene_overlap_arr[$i]}
    qmax_overlap=${qmax_overlap_arr[$i]}
    bp_win=${bp_win_arr[$i]}
    min_ct=${min_ct_arr[$i]}
    
    name=${cohort}_${bin_size}_${min_anchor}_${max_gap}_${max_overlap}_${min_gene_gap}_${max_gene_overlap}_${qmax_overlap}_${bp_win}_${min_ct}_${filt}
    fd_output_path=/nas/longleaf/home/jclin/jwanglab_users_jclin/data/10_fusion_detection/07_new_fd/${name}

    echo $name
    log_path=$fd_output_path/logs

    mkdir -p $fd_output_path
    mkdir -p $log_path

    deps="--dependency=afterok"
    {
    read # skip the first line ie header
    while IFS=$'\t' read -r -a array; do
        seq_id=${array[$seq_id_col]}
        echo "fusilli processing $seq_id..."
       id=$(sbatch -n1 -N1 -c1 --mem=16G --time=11-0 --job-name=01_fusilli_${seq_id} --out=${log_path}/01_fusilli_${seq_id}.log --wrap="python3 /proj/jwanglab/users/jclin/nanopore-dx/10_fusion_detection/fusilli/src/fusilli/fusilli.py -p ${input_path}/${seq_id}.paf -o ${fd_output_path} -r -nfm")
        id=${id##* }
        deps="$deps:$id"
    done
    } < $metadata_path
    echo $deps

    sbatch ${deps} -n1 -N1 -c1 --mem=16G --time=11-0 --job-name=02_${name}_aggregate_results --out=${log_path}/02_${name}_aggregate_results.log --wrap="python3 /nas/longleaf/home/jclin/jwanglab_users_jclin/nanopore-dx/10_fusion_detection/fusilli/src/fusilli/aggregate_fusilli.py ${fd_output_path} ${fusion_master_path} ${fusion_std_path} ${name} ${metadata_path} ${metadata_fusion} ${bin_size} ${min_anchor} ${max_gap} ${max_overlap} ${min_gene_gap} ${max_gene_overlap} ${qmax_overlap} ${bp_win} ${min_ct} ${cohort} ${filt} ${comment} ${gu_fusion_std_path} ${gu_pax5_status_std_path} ${wlist_path2} > ${fd_output_path}/${name}.log"
done

