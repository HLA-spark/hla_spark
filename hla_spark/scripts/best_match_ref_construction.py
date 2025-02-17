#!/root/anaconda3/bin/python
# -*- coding: UTF-8 -*-

## CreateTime: 2023-11-30  14:01


import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

def high_aligment_score_bam(time, dir_fastq, bam, samtools, as_threshold, ncpu):
    key_sam = f"{dir_fastq}/{time}.sam"
    filter_sam = f"{dir_fastq}/{time}.AS.sam"
    filter_bam = f"{dir_fastq}/{time}.AS.bam"
    os.system(
        '{samtools} view -@ {ncpu} -h {bam} > {key_sam}'.format(bam=bam, key_sam=key_sam, samtools=samtools, ncpu=ncpu))
    high_aligment_score_sam(samtools, as_threshold, key_sam, filter_sam, filter_bam, ncpu)
    return filter_bam


def high_aligment_score_sam(samtools, as_threshold, key_sam, filter_sam, filter_bam, ncpu):
    f_in = open(key_sam, 'r')
    f_out = open(filter_sam, 'w')
    lines = f_in.readlines()
    for line in lines:
        line = line.strip()
        if '@' in line:
            f_out.write(line + '\n')
        else:
            tokens = line.split('\t')
            mapq = int(tokens[4])
            tokens = line.split('AS:i:')
            as_i = int(tokens[1].split('\t')[0])

            if (as_i > as_threshold) and (mapq != 255):
                f_out.write(line + '\n')

    f_in.close()
    f_out.close()
    os.system(f'{samtools} view -bS {filter_sam} -@ {ncpu} | samtools sort -@ {ncpu} -o {filter_bam}')
    os.system(f'{samtools} index {filter_bam}')
    os.system(f'rm {key_sam}')
    os.system(f'rm {filter_sam}')
    return filter_bam


def coverage_cal_exon23(hla_current, bed, filter_bam, dir_fastq, df_ref_list,bedtools):
    filter_depth = f"{dir_fastq}/{hla_current}.AS.depth"
    os.system(f"{bedtools} coverage -a {bed} -b {filter_bam}>{filter_depth}")
    df = pd.read_csv(filter_depth, sep='\t', header=None)
    df.columns = ['chr', 'start', 'end', 'name', 'reads', 'bases', 'size', 'coverage']
    df = df[df['coverage'] != 0]
    df['exon'] = df['name'].apply(lambda x: x[0:4])
    df_stat = pd.DataFrame(columns=['cov_exon23', 'depth_exon23', 'cov_exon', 'depth_exon'])
    chr_set = set(df['chr'])
    df_exon23 = df[(df['name'] == 'exon2') | (df['name'] == 'exon3')]
    df_exon = df[df['exon'] == 'exon']
    for i in chr_set:
        cov_exon23 = np.sum(df_exon23[df_exon23['chr'] == i]['coverage'])
        depth_exon23 = np.sum(df_exon23[df_exon23['chr'] == i]['reads'])
        cov_exon = np.sum(df_exon[df_exon['chr'] == i]['coverage'])
        depth_exon = np.sum(df_exon[df_exon['chr'] == i]['reads'])
        df_stat.loc[i] = [cov_exon23, depth_exon23, cov_exon, depth_exon]
    df_max = df_stat[df_stat['cov_exon23'] != 0]
    if len(df_max.index) == 0:  ##第一种情况：没有找到任何相似单倍型，用0标注
        return 0, None, 0, None
    df_max = df_max.join(df_ref_list, how='left')
    df_max = df_max.sort_values(by=['cov_exon23', 'depth_exon23', 'cov_exon', 'depth_exon'],
                                ascending=False)
    max_cov_exon23 = df_max.iloc[0, 0]
    max_num = len(df_max.index)
    df_max.to_csv(dir_fastq + '/countlist.' + hla_current + '.xls', sep='\t')
    return max_num, df_max.index.tolist(), max_cov_exon23, df_max


def coverage_cal_exon(hla_current, bed, filter_bam, dir_fastq, df_ref_list,bedtools):
    filter_depth = f"{dir_fastq}/{hla_current}.AS.depth"
    os.system(f"{bedtools} coverage -a {bed} -b {filter_bam}>{filter_depth}")
    df = pd.read_csv(filter_depth, sep='\t', header=None)
    df.columns = ['chr', 'start', 'end', 'name', 'reads', 'bases', 'size', 'coverage']
    df = df[df['coverage'] != 0]
    df['exon'] = df['name'].apply(lambda x: x[0:4])
    df_stat = pd.DataFrame(columns=['cov_exon23', 'depth_exon23', 'cov_exon', 'depth_exon'])
    chr_set = set(df['chr'])
    df_exon23 = df[(df['name'] == 'exon2') | (df['name'] == 'exon3')]
    df_exon = df[df['exon'] == 'exon']
    for i in chr_set:
        cov_exon23 = np.sum(df_exon23[df_exon23['chr'] == i]['coverage'])
        depth_exon23 = np.sum(df_exon23[df_exon23['chr'] == i]['reads'])
        cov_exon = np.sum(df_exon[df_exon['chr'] == i]['coverage'])
        depth_exon = np.sum(df_exon[df_exon['chr'] == i]['reads'])
        df_stat.loc[i] = [cov_exon23, depth_exon23, cov_exon, depth_exon]
    df_max = df_stat[df_stat['cov_exon'] != 0]
    if len(df_max.index) == 0:  ##第一种情况：没有找到任何相似单倍型，用0标注
        return 0, None, 0, None
    df_max = df_max.join(df_ref_list, how='left')
    df_max = df_max.sort_values(by=['cov_exon', 'depth_exon', 'cov_exon23', 'depth_exon23'],
                                ascending=False)
    max_cov_exon = df_max.iloc[0, 2]
    # 如果多个则进入有剔除的迭代筛选
    max_num = len(df_max.index)
    df_max.to_csv(dir_fastq + '/countlist.' + hla_current + '.xls', sep='\t')
    return max_num, df_max.index.tolist(), max_cov_exon, df_max


def get_unmapped(sam_file, sam_unmapped, as_threshold, bam_unmapped, samtools, ncpu):
    f_in = open(sam_file, 'r')
    f_out = open(sam_unmapped, 'w')
    lines = f_in.readlines()
    for line in lines:
        line = line.strip()
        if '@' in line:
            f_out.write(line + '\n')
        else:
            tokens = line.split('AS:i:')
            as_i = int(tokens[1].split('\t')[0])
            if (as_i > 145):
                continue
            else:
                f_out.write(line + '\n')
    f_in.close()
    f_out.close()
    # 转化为sorted.bam
    os.system(f'{samtools} view -bS {sam_unmapped} -@ {ncpu} | {samtools} sort -@ {ncpu} -o {bam_unmapped}')
    os.system(f'{samtools} index {bam_unmapped}')
    return


def unmapped_extract(df_ref_list, dir_fastq, hla_current, bam, bedtools, bwa, sample, ncpu, as_threshold, samtools):
    fasta_file = f"{dir_fastq}/{hla_current}.fasta"
    fq1 = f"{dir_fastq}/{hla_current}.r1.fq"
    fq2 = f"{dir_fastq}/{hla_current}.r2.fq"
    sam_file = f"{dir_fastq}/{hla_current}.align.sam"
    sam_unmapped = f"{dir_fastq}/{hla_current}.umapped.sam"
    bam_unmapped = f"{dir_fastq}/{hla_current}.unmapped.bam"
    fq1_un = f"{dir_fastq}/{hla_current}.unmapped.r1.fq"
    fq2_un = f"{dir_fastq}/{hla_current}.unmapped.r2.fq"
    # 提取current hla的fasta
    write_fasta([hla_current], df_ref_list, fasta_file,bwa)
    # bam转化为fastq
    os.system(f'{bedtools} bamtofastq -i {bam} -fq {fq1} -fq2 {fq2}')
    # 将bam比对到hla
    os.system(
        f'{bwa} mem -M -R "@RG\\tID:{sample}\\tSM:{sample}\\tLB:genome\\tPL:ILLUMINA" -t {ncpu} {fasta_file} {fq1} {fq2} -o {sam_file}')
    ##提取unmapped reads
    get_unmapped(sam_file, sam_unmapped, as_threshold, bam_unmapped, samtools, ncpu)
    os.system(f'{bedtools} bamtofastq -i {bam_unmapped} -fq {fq1_un} -fq2 {fq2_un}')

    os.system(f'rm {sam_file}')
    os.system(f'rm {sam_unmapped}')
    os.system(f'rm {fasta_file}*')
    os.system(f'rm {fq1}')
    os.system(f'rm {fq2}')
    os.system(f'rm {bam_unmapped}*')

    return fq1_un, fq2_un


def align_sub(hla_current, hla_others, df_ref_list, dir_fastq, sample, ncpu, as_threshold, samtools, bwa, fq1_un,
              fq2_un):  # 1.提取current hla的fasta，将bam比对到hla，截取未比及比对质量小于130的reads，生成sort bam；
    fasta_others = f"{dir_fastq}/{hla_current}.others.fasta"
    sam_others = f"{dir_fastq}/{hla_current}.other.sam"
    sam_others_filter = f"{dir_fastq}/{hla_current}.other.filter.sam"
    bam_others = f"{dir_fastq}/{hla_current}.other.bam"

    ##比对到hla_others
    write_fasta(hla_others, df_ref_list, fasta_others,bwa)
    os.system(
        f'{bwa} mem -M -R "@RG\\tID:{sample}\\tSM:{sample}\\tLB:genome\\tPL:ILLUMINA" -t {ncpu} {fasta_others} {fq1_un} {fq2_un} -o {sam_others}')

    bam_others = high_aligment_score_sam(samtools, as_threshold, sam_others, sam_others_filter, bam_others, ncpu)

    os.system(f'rm {fasta_others}*')
    return bam_others


def ref_iterator(hla_types_total, df_ref_list, dir_fastq, bam, bed, sample, ncpu, as_threshold, samtools, bedtools, bwa,
                 is_hla_gene):
    hla_types_selected = [hla_types_total[0]]
    hla_others = hla_types_total[1:]
    hla_current = hla_types_selected[0]
    fq1_un, fq2_un = unmapped_extract(df_ref_list, dir_fastq, hla_current, bam, bedtools, bwa, sample, ncpu,
                                      as_threshold, samtools)
    j = '2nd'
    current_bam = align_sub(j, hla_others, df_ref_list, dir_fastq, sample, ncpu, as_threshold, samtools, bwa, fq1_un,
                            fq2_un)
    if is_hla_gene == 1:
        max_num, hla_types, max_cov, df_max = coverage_cal_exon23(j, bed, current_bam, dir_fastq, df_ref_list,bedtools)
    else:
        max_num, hla_types, max_cov, df_max = coverage_cal_exon(j, bed, current_bam, dir_fastq, df_ref_list,bedtools)
    if max_cov > 0:  ##第一种情况：没有找到任何相似单倍型，用0标注
        hla_types_selected.append(hla_types[0])
    os.system(f'rm {fq1_un}')
    os.system(f'rm {fq2_un}')
    return hla_types_selected


def write_fasta(hla_types, df_ref_list, fasta_file,bwa):
    out = open(fasta_file, 'w')
    for h in hla_types:
        out.write('>' + h + ' ' + df_ref_list.loc[h, 'name_full'] + ' ' + str(
            df_ref_list.loc[h, 'length']) + '\n')
        out.write(df_ref_list.loc[h, 'sequence'] + '\n')
    out.close()
    os.system(f'{bwa} index {fasta_file}')


def ref_construction(hla_types, df_ref_list, ref, bestmatch, max_cov, sample,bwa):
    write_fasta(hla_types, df_ref_list, ref,bwa)
    out = open(bestmatch, 'w')
    out.write('\t'.join([sample, ';'.join(hla_types), str(max_cov)]))
    out.close()

def run(dir_fastq,sample,bed,bam,ref,ncpu,is_hla_gene,config,data_path):  #dir, sample, gene_bed, gene_sort_bam, gene_fasta, ncpu, is_hla_gene
    reflist = data_path+config.get('paths', 'HLA_alleles_pool')
    samtools = config.get('soft', 'Samtools')
    bedtools = config.get('soft', 'Bedtools')
    bwa = config.get('soft', 'BWA')
    as_threshold = int(config.get('paras', 'alignment_score_threshold'))

    df_ref_list = pd.read_csv(reflist, sep='\t', index_col=0)
    bestmatch = f'{dir_fastq}/bestmatch'
    filter_bam = high_aligment_score_bam("1st", dir_fastq, bam, samtools, as_threshold, ncpu)
    if is_hla_gene == 1:
        max_num, hla_types, max_cov, df_max = coverage_cal_exon23('1st', bed, filter_bam, dir_fastq, df_ref_list,bedtools)
    else:
        max_num, hla_types, max_cov, df_max = coverage_cal_exon('1st', bed, filter_bam, dir_fastq, df_ref_list,bedtools)

    if max_cov == 0:
        out = open(bestmatch, 'w')
        out.write(sample + '\t\n')
        out.close()
    elif max_num == 1:
        ref_construction(hla_types, df_ref_list, ref, bestmatch, max_cov, sample,bwa)
    else:
        hla_types = ref_iterator(hla_types, df_ref_list, dir_fastq, bam, bed, sample, ncpu, as_threshold, samtools,
                                 bedtools, bwa, is_hla_gene)
        ref_construction(hla_types, df_ref_list, ref, bestmatch, max_cov, sample,bwa)