#!/usr/bin/python
# -*- coding: UTF-8 -*-
## FileName  : HLA_spark.py


import os
current_dir=os.path.dirname(os.path.abspath(__file__))
data_path=os.path.join(current_dir,"data")+'/'
configfile=os.path.join(current_dir,'config','config_HLA_spark.ini')

import configparser
import pandas as pd
from multiprocessing import Pool
import warnings
warnings.filterwarnings("ignore")
import subprocess

from hla_spark.scripts import filter_by_alignment_score as filter_al, best_match_ref_construction as best_m, \
    genomic_coordinate_transformation as coord_stand
from hla_spark.scripts.common import *
from hla_spark.scripts.argument_parser import *


sample=None
fq1=None
fq2=None
outdir=None
ncpu=4
config=None
steps_list=None

dir_hla_reads=None
hla_bam_sort=None
dir_snp_calling=None
df_genelist=None
dir_snp_report=None
sub_log=None

def arg_init1(sample_name, fastq1, fastq2, output_dir, n,step_list):
    global sample
    global fq1
    global fq2
    global outdir
    global ncpu
    global config
    global sub_log
    global steps_list
    global configfile

    sample = sample_name
    fq1 = fastq1
    fq2=fastq2
    ncpu=n
    config = configparser.ConfigParser()
    config.read(configfile)
    outdir = mkdir(output_dir)
    sub_log_dir = mkdir(f"{outdir}/{sample_name}/Log")
    sub_log = backup_std_err(sample_name, sub_log_dir)
    steps_list = step_list.split(",")

def arg_init2():
    global dir_hla_reads
    global hla_bam_sort
    global dir_snp_calling
    global df_genelist
    global dir_snp_report
    dir_hla_reads = mkdir(f"{outdir}/{sample}/hla_reads")
    hla_bam_sort = f"{dir_hla_reads}/{sample}_hla_typing.sort.bam"

    dir_snp_calling = mkdir(f"{outdir}/{sample}/snp_calling")

    df_genelist = pd.read_csv(data_path+config.get("paths", "HLA_genes"), sep='\t', index_col=0)

    dir_snp_report = mkdir(f"{outdir}/{sample}/snp_report")


def HLA_reads_extraction_module():
    cmd = []
    cmd.append('\n------------------------------Mapping reads to the HLA reference sequence database----------------------------')
    print('------------------------------Mapping reads to the HLA reference sequence database----------------------------')
    fq1_trim = f"{dir_hla_reads}/{sample}.qc.r1.fq.gz"
    fq2_trim = f"{dir_hla_reads}/{sample}.qc.r2.fq.gz"
    start_time = startTime()
    if config.get("paras", 'qc') == 'yes':
        print('Fastp analysis in process ...')
        fq1_tmp = f"{dir_hla_reads}/{sample}.qc.r1.fq"
        fq2_tmp = f"{dir_hla_reads}/{sample}.qc.r2.fq"
        json_trim = f"{dir_hla_reads}/{sample}_fastp_report.json"
        html_trim = f"{dir_hla_reads}/{sample}_fastp_report.html"
        command="{fastp} --in1 {fq1} --in2 {fq2} --out1 {fq1_trim} --out2 {fq2_trim} -j {json_trim} -h {html_trim}  -3 -W 4 -M 15 -q 15 -l 36 -w {ncpu}".format(fastp=config.get("soft", "Fastp"), fq1=fq1, fq2=fq2, fq1_trim=fq1_tmp,
                 fq2_trim=fq2_tmp, json_trim=json_trim, html_trim=html_trim, ncpu=ncpu)
        subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
        rmfile(fq1_trim)
        rmfile(fq2_trim)
        command='{pigz} --best {fq1_new}'.format(pigz=config.get('soft', 'Pigz'), fq1_new=fq1_tmp)
        subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
        command='{pigz} --best {fq2_new}'.format(pigz=config.get('soft', 'Pigz'), fq2_new=fq2_tmp)
        subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
    else:
        command='ln -s {fq1} {fq1_trim}'.format(fq1=fq1, fq1_trim=fq1_trim)
        subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
        command='ln -s {fq2} {fq2_trim}'.format(fq2=fq2, fq2_trim=fq2_trim)
        subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)

    hla_bam = f"{dir_hla_reads}/{sample}_hla_typing.bam"
    print('BWA mem analysis in process ...')
    command='{bwa} mem -M -R "@RG\\tID:{sample}\\tSM:{sample}\\tLB:genome\\tPL:ILLUMINA" -t {ncpu} {hla_ref} {fq1_new} {fq2_new}|{samtools} view -Sb -bF 12 -@ {ncpu} -o - >{hla_bam}'.format(bwa=config.get("soft", "BWA"), ncpu=ncpu, hla_ref=data_path+config.get('paths', 'HLA_ref'),
                fq1_new=fq1_trim, fq2_new=fq2_trim, sample=sample,
                samtools=config.get("soft", "Samtools"), hla_bam=hla_bam)
    subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
    command='{samtools} sort -O bam -@ {ncpu} -o {hla_bam_sort} {hla_bam}'.format(samtools=config.get("soft", "Samtools"), ncpu=ncpu, hla_bam_sort=hla_bam_sort,hla_bam=hla_bam)
    subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
    command='{samtools} index {hla_bam_sort}'.format(samtools=config.get("soft", "Samtools"),
                                                        hla_bam_sort=hla_bam_sort)
    subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
    rmfile(hla_bam)
    spend_time = spendTime(start_time)
    cmd.append(f'\nJob of \"Mapping reads to the HLA reference sequence database\" is completed. Total Spend time is {spend_time}\n')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()
    return

def getlist():
    f = open(data_path+config.get('paths', 'HLA_exon23_list'), 'r')
    exon23 = []
    lines = f.readlines()
    for line in lines:
        exon23.append(line.strip())
    f.close()
    return exon23


def Best_matched_alleles_module(gene):
    exon23 = getlist()
    cmd = []
    cmd.append(
        f'\n------------------------------Best-matching HLA allele selection for HLA gene: {gene}----------------------------')
    print(
        f'------------------------------Best-matching HLA allele selection for HLA gene: {gene}----------------------------')
    dir_gene = mkdir(f"{dir_snp_calling}/{gene}/")
    dir = mkdir(f"{dir_gene}/00.fastq")
    is_hla_gene = 0
    if gene in exon23:
        is_hla_gene = 1
    gene_bed = data_path+config.get('paths', 'HLA_bed') + '/' + gene + '.bed'
    gene_fasta = f"{dir}/{gene}_bestm.fasta"
    gene_bam = f"{dir}/{sample}.{gene}.bam"
    gene_sort_bam = f"{dir}/{sample}.{gene}.sort.bam"
    fq1 = f"{dir}/{sample}.{gene}.r1.fq"
    fq2 = f"{dir}/{sample}.{gene}.r2.fq"
    cmd.append('\nRead binning in process ...')
    print('Read binning in process ...')
    command='{samtools} view -L {gene_bed} -Sb {hla_bam} > {gene_bam}'.format(samtools=config.get("soft", "Samtools"), gene_bed=gene_bed,
                      hla_bam=hla_bam_sort,gene_bam=gene_bam)
    subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
    command='{samtools} sort -n {gene_bam} -o {gene_sort_bam}'.format(samtools=config.get("soft", "Samtools"), gene_bam=gene_bam,
                      gene_sort_bam=gene_sort_bam)
    subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
    command='{bedtools} bamtofastq -i {gene_sort_bam} -fq {fq1_pair} -fq2 {fq2_pair}'.format(bedtools=config.get("soft", "Bedtools"), gene_sort_bam=gene_sort_bam, fq1_pair=fq1,
                       fq2_pair=fq2)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    #subprocess.call(command,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    rmfile(f'{fq1}.gz')
    rmfile(f'{fq2}.gz')
    command='{pigz} --best {fq1}'.format(pigz=config.get('soft', 'Pigz'), fq1=fq1)
    subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
    command='{pigz} --best {fq2}'.format(pigz=config.get('soft', 'Pigz'), fq2=fq2)
    subprocess.call(command,shell=True,stdout=sub_log,stderr=sub_log)
    rmfile(gene_bam)
    cmd.append('\nSelecting best-match alleles in process ...')
    print('Selecting best-match alleles in process ...')
    best_m.run(dir, sample, gene_bed, gene_sort_bam, gene_fasta, ncpu, is_hla_gene,config,data_path)
    cmd.append(f'\nJob of \"Best-matching HLA allele selection\" for HLA gene: {gene} is completed.\n')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()
    return


def Align_module(gene):
    cmd = []
    cmd.append(
        f"\n-------------------------------Alignment to best-matching HLA allele reference sequence: {gene}-------------------------------")
    print(f"-------------------------------Alignment to best-matching HLA allele reference sequence: {gene}-------------------------------")
    dir_gene = mkdir(f"{dir_snp_calling}/{gene}/")
    dir_fastq = mkdir(f"{dir_gene}/00.fastq")
    dir_bam = mkdir(f"{dir_gene}/01.bam")

    fq1 = f"{dir_fastq}/{sample}.{gene}.r1.fq.gz"
    fq2 = f"{dir_fastq}/{sample}.{gene}.r2.fq.gz"
    aln_sam = f"{dir_bam}/{sample}.{gene}.aln.sam"
    qc_sam = f"{dir_bam}/{sample}.{gene}.qc.sam"
    qc_bam = f"{dir_bam}/{sample}.{gene}.qc.bam"
    aln_sort_bam = f"{dir_bam}/{sample}.{gene}.aln.sort.bam"
    gene_fasta = f"{dir_fastq}/{gene}_bestm.fasta"
    command=config.get("soft", "BWA") + ' index ' + gene_fasta
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    command=config.get("soft", "Samtools") + ' faidx ' + gene_fasta
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    cmd.append('\nBWA mem analysis in process ...')
    print('BWA mem analysis in process ...')
    command='{bwa} mem -M -R "@RG\\tID:{sample}\\tSM:{sample}\\tLB:genome\\tPL:ILLUMINA" -t {ncpu} {ref} {fq1} {fq2} -o {aln_sam}'.format(
        bwa=config.get("soft", "BWA"), ncpu=ncpu, ref=gene_fasta,fq1=fq1, fq2=fq2, sample=sample, aln_sam=aln_sam)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    print('Extracting reads with high alignment scores ...')
    filter_al.run(aln_sam, config.get('paras', 'alignment_score_threshold'), qc_sam)
    command="{samtools} view -bh -F 4 -F 256 {qc_sam} -o {qc_bam}".format(samtools=config.get("soft", "Samtools"), qc_sam=qc_sam, qc_bam=qc_bam)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    command= "{samtools} sort -@ {ncpu} -o {aln_sort_bam} {qc_bam}".format(samtools=config.get("soft", "samtools"), ncpu=ncpu, qc_bam=qc_bam, aln_sort_bam=aln_sort_bam)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    command="{samtools} index {aln_sort_bam}".format(samtools=config.get("soft", "Samtools"), aln_sort_bam=aln_sort_bam)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    rmfile(aln_sam)
    rmfile(qc_sam)
    rmfile(qc_bam)
    cmd.append(f'\nJob of \"Alignment to best-matching HLA allele reference sequence\" for HLA gene: {gene} is completed.\n')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()
    return


def Rmdup_module(gene):
    cmd = []
    cmd.append(
        f"\n-------------------------------Removing duplication reads from the BAM file: {gene}-------------------------------")
    print(f"-------------------------------Removing duplication reads from the BAM file: {gene}-------------------------------")
    dir_gene = mkdir(f"{dir_snp_calling}/{gene}/")
    dir_bam = mkdir(f"{dir_gene}/01.bam")
    aln_sort_bam = f"{dir_bam}/{sample}.{gene}.aln.sort.bam"

    rmdup_bam = f"{dir_bam}/{sample}.{gene}.rmdup.bam"  #
    rmdup_log = f"{dir_bam}/{sample}.{gene}.rmdup.log"  #
    rmdup_sorted_bam = f"{dir_bam}/{sample}.{gene}.rmdup.sort.bam"
    cmd.append('\nMarking duplicates by picard ...')
    print('Marking duplicates by picard ...')
    command="{picard} MarkDuplicates  I={aln_sort_bam} O={rmdup_bam} M={rmdup_log}".format \
            (picard=config.get("soft", "Picard"), aln_sort_bam=aln_sort_bam, rmdup_bam=rmdup_bam,rmdup_log=rmdup_log)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    print('Sorting and indexing by samtools ...')
    command="{samtools} sort -@ {ncpu} -o {rmdup_sorted_bam} {rmdup_bam}".format \
            (samtools=config.get("soft", "samtools"), ncpu=ncpu, rmdup_bam=rmdup_bam,rmdup_sorted_bam=rmdup_sorted_bam)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    command="{samtools} index {rmdup_sorted_bam}".format \
            (samtools=config.get("soft", "samtools"), rmdup_sorted_bam=rmdup_sorted_bam)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    if os.path.exists(rmdup_sorted_bam):
        rmfile(aln_sort_bam)
        rmfile(rmdup_bam)
    cmd.append(
        f'\nJob of \"Removing duplication reads from the BAM file\" for HLA gene: {gene} is completed.\n')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()
    return


def Mpileup_module(gene):
    cmd = []
    cmd.append(
        f"\n-------------------------------Calculating genotype likelihoods by BCFtools mpileup: {gene}-------------------------------")
    print(
        f"-------------------------------Calculating genotype likelihoods by BCFtools mpileup: {gene}-------------------------------")
    dir_gene = mkdir(f"{dir_snp_calling}/{gene}/")
    dir_fastq = mkdir(f"{dir_gene}/00.fastq")
    dir_bam = mkdir(f"{dir_gene}/01.bam")
    dir_mpileup = mkdir(f"{dir_gene}/02.mpileup")
    gene_fasta = f"{dir_fastq}/{gene}_bestm.fasta"
    rmdup_sorted_bam = f"{dir_bam}/{sample}.{gene}.rmdup.sort.bam"
    mpileup = f"{dir_mpileup}/{sample}.{gene}.mpileup.vcf"
    command="{bcftools} mpileup -a FORMAT/DP,FORMAT/AD {rmdup_sorted_bam} -A --fasta-ref {ref} > {mpileup}".format \
            (bcftools=config.get('soft', 'Bcftools'), ref=gene_fasta, rmdup_sorted_bam=rmdup_sorted_bam,mpileup=mpileup)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    cmd.append(
        f'\nJob of \"Calculating genotype likelihoods by BCFtools mpileup\" for HLA gene: {gene} is completed.\n')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()
    return


def VCFConversion_module(gene):
    cmd = []
    cmd.append(
        f"\n-------------------------------Genomic coordinate standardization based on GRCh38.p14: {gene}-------------------------------")
    print(
        f"-------------------------------Genomic coordinate standardization based on GRCh38.p14: {gene}-------------------------------")
    dir_gene = mkdir(f"{dir_snp_calling}/{gene}/")
    dir_mpileup = mkdir(f"{dir_gene}/02.mpileup")
    mpileup = f"{dir_mpileup}/{sample}.{gene}.mpileup.vcf"
    mut_hg38 = f"{dir_mpileup}/{sample}.{gene}.genotype.hg38.vcf"
    hla_matrix = data_path+config.get('paths', 'HLA_matrix_path') + '/' + gene + '.xls'
    chrom = df_genelist.loc[gene, 'chr']
    coord_stand.run(mpileup, hla_matrix, mut_hg38, chrom, sample)
    cmd.append(
        f'\nJob of \"Genomic coordinate standardization based on GRCh38.p14\" for HLA gene: {gene} is completed.\n')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()
    return


def Genotype_report_module():
    # vcf combine by gene
    cmd = []
    headers = []
    headers.append('##fileformat=VCFv4.2')
    headers.append('##FILTER=<ID=PASS,Description="All filters passed">')
    headers.append('##contig=<ID=chr6,length=170805979>')
    headers.append('##contig=<ID=HSCHR6_MHC_QBL_CTG1,length=4606388>')
    headers.append('##contig=<ID=HSCHR6_MHC_SSTO_CTG1,length=4929269>')
    headers.append('##ALT=<ID=*,Description="Represents allele(s) other than observed.">')
    headers.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">')
    headers.append(
        '##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">')
    headers.append('##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">')
    headers.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    headers.append('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">')
    headers.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">')
    headers.append('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">')
    headers.append(
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample]) + '\n')
    header_file = dir_snp_report + '/header'
    dat_file = dir_snp_report + '/dat'
    vcf_file = dir_snp_report + '/' + sample + '.genotype.unsort.vcf'
    vcf_sort = dir_snp_report + '/' + sample + '.genotype.vcf'
    f = open(header_file, 'w')
    f.write('\n'.join(headers))
    f.close()
    cmd.append(
        '\n------------------------------Merging genotype VCF files by genes----------------------------')
    print(
        '------------------------------Merging genotype VCF files by genes----------------------------')
    start_time = startTime()
    command='cat ' + dir_snp_calling + '/*/02.mpileup/*genotype.hg38.vcf>' + dat_file
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    command='cat ' + header_file + ' ' + dat_file + '>' + vcf_file
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    rmfile(dat_file)
    command='{picard} SortVcf I={vcf} O={vcf_sort}'.format(picard=config.get('soft', 'Picard'), vcf=vcf_file,vcf_sort=vcf_sort)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    gz_vcf = vcf_sort + '.gz'
    command=config.get('soft', 'Bgzip') + ' -c -f ' + vcf_sort + ' > ' + gz_vcf
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    command='{bcftools} index {gz_vcf}'.format(bcftools=config.get('soft', 'Bcftools'), gz_vcf=gz_vcf)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    rmfile(vcf_file)
    spend_time = spendTime(start_time)
    cmd.append(
        f'\nJob of \"Merging genotype VCF files by genes\" is completed. Total Spend time is {spend_time}\n')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()
    return


def SNP_calling_module():
    cmd = []
    cmd.append(
        '\n------------------------------Calling SNPs in HLA genes by BCFtools call----------------------------')
    print(
        '------------------------------Calling SNPs in HLA genes by BCFtools call----------------------------')
    start_time = startTime()
    genotype_sort = dir_snp_report + '/' + sample + '.genotype.vcf.gz'
    snpfile = dir_snp_report + '/' + sample + '.snp.vcf'
    command='{bcftools} call -O v --threads {ncpu} -m --ploidy 2 -p 0.05 -o {snpfile} {genotype_sort}'.format(
        bcftools=config.get('soft', 'Bcftools'), ncpu=ncpu, snpfile=snpfile, genotype_sort=genotype_sort)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    gz_vcf = snpfile + '.gz'
    command=config.get('soft', 'Bgzip') + ' -c -f ' + snpfile + ' > ' + gz_vcf
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    command='{bcftools} index {gz_vcf}'.format(bcftools=config.get('soft', 'Bcftools'), gz_vcf=gz_vcf)
    subprocess.call(command, shell=True, stdout=sub_log, stderr=sub_log)
    spend_time = spendTime(start_time)
    cmd.append(
        f'\nJob of \"Calling SNPs in HLA genes by BCFtools call\" is completed. Total Spend time is {spend_time}\n')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()
    return

def gene_jobs(gene):
    cmd = []

    cmd.append(
        f'\n------------------------------HLA gene: {gene} in analysis----------------------------')
    print(
        f'------------------------------HLA gene: {gene} in analysis----------------------------')
    start_time = startTime()
    if "2" in steps_list:
        Best_matched_alleles_module(gene)
    gene_fasta = f"{dir_snp_calling}/{gene}/00.fastq/{gene}_bestm.fasta"
    print(gene_fasta)
    if os.path.exists(gene_fasta):
        if "3" in steps_list:
            Align_module(gene)
        if "4" in steps_list:
            Rmdup_module(gene)
        if "5" in steps_list:
            Mpileup_module(gene)
        if "6" in steps_list:
            VCFConversion_module(gene)

    spend_time = spendTime(start_time)
    cmd.append(
        f'\nJobs of HLA gene: {gene} is completed. Total Spend time is {spend_time}')
    sub_log.write('\n'.join(cmd))
    sub_log.flush()

def run_hla_spark(sample_name, fastq1, fastq2, output_dir="hla_spark_results/", n_core=4,step_list="1,2,3,4,5,6,7,8",hla_genes='all'):
    arg_init1(sample_name,fastq1,fastq2,output_dir,n_core,step_list)
    arg_init2()
    try:
        if hla_genes=='all':
            hla_genes = list(df_genelist.index)
        else:
            hla_genes=hla_genes.split(',')
        print('Sample '+sample+' analysis is started......')
        if "1" in steps_list:
            HLA_reads_extraction_module()
        if ("2" in steps_list) or ("3" in steps_list) or ("4" in steps_list) or ("5" in steps_list) or (
                "6" in steps_list):
            pool=Pool(processes=ncpu)
            res=pool.map(gene_jobs,hla_genes)
            pool.close()
            pool.join()
        if "7" in steps_list:
            Genotype_report_module()
        if "8" in steps_list:
            SNP_calling_module()
    finally:
        sub_log.flush()
        sub_log.close()

def hla_spark():
    arguments = argument_parser()
    run_hla_spark(**arguments)

if __name__ == '__main__':
    arguments = argument_parser()
    run_hla_spark(**arguments)
