# HLA-spark: An accurate method for HLA SNP calling facilitates prioritization of genetic factors on occult hepatitis B virus infection at SNP-level

***

**HLA-spark** is a method for identifying SNPs in human leukocyte antigen (HLA) genes from next-generation sequencing (NGS) data. It provides reliable genotype estimates by using the best-matching HLA allele(s) as the reference and transforming the HLA allele-based genotype file to GRCh38.p14 genomic coordinates. **HLA-spark** addresses the challenges of directly calling HLA mutations from NGS data, which can then be utilized in phenotype-association analyses.

![graphical_abstract](https://github.com/user-attachments/assets/635530ed-bfb3-404c-a1e3-ba079ef1d3aa)

[**hla\_spark**](https://github.com/HLA-spark/hla_spark) is an HLA SNP caller for analysing a sample data. 

[**hla\_spark\_merge**](https://github.com/HLA-spark/hla_spark_merge) is used for merging and annotating SNP results of multiple samples to facilitate the further phenotype-association analysis.

## **hla\_spark** Installation

***

Download the package from Github:

    git clone https://github.com/HLA-spark/hla_spark.git
    cd hla_spark/hla_spark

Download data.tar.gz from Github in the path of **hla_spark/hla_spark/**, and extract the data files 

    wget -c https://github.com/HLA-spark/raw/refs/heads/main/hla_spark/data.tar.gz?download=
    tar -zxvf data.tar.gz

Install Anaconda if you do not have it.

Create the conda environment with environment.yml, and then activate it as follows in terminal.

    cd ../
    conda env create -f environment.yml -n hla_spark
    conda activate hla_spark

Get the path of conda on Linux

```
whereis conda

```

>     conda: /root/anaconda3/bin/conda

Install HLA-spark

    /root/anaconda3/envs/hla_spark/bin/pip install .

## **hla\_spark** usage&#x20;

***

### Running **hla\_spark**

Check the command line options

    hla_spark --help

### Arguments for **hla\_spark**

Required arguments

1.  \-s/--sample\_name: used for naming the analysis direcotry. A subdirectory named the input sample name will be auto-created under the output-dir.
2.  \-f1/--fastq1: the file path of r1.fq(.gz) data
3.  \-f2/--fastq2: the file path of r2.fq(.gz) data

Optional arguments

1.  \-o/--output-dir, default="hla\_spark\_results/": Directory path for the output. You can put the analysis results of multiple samples under a same directory. A subdirectory named the input sample name will be auto-created under the output-dir.
2.  \-n/--n-core, default=4: Number of threads to use in the analysis
3.  \-l/--step-list, default="1,2,3,4,5,6,7,8":&#x20;

    **hla\_spark** contains eight analysis modules. Analysis is permitted at selected steps considering the unexcepted interruption.

    **1 :** Mapping reads to the HLA reference sequence database&#x20;

    **2 :** Best-matching HLA allele selection for HLA gene&#x20;

    **3 :** Alignment to best-matching HLA allele reference sequence&#x20;

    **4 :** Removing duplication reads from the BAM file&#x20;

    **5 :** Calculating genotype likelihoods&#x20;

    **6 :** Genomic coordinate standardization based on GRCh38.p14&#x20;

    **7 :** Merging genotype VCF files by genes&#x20;

    **8 :** Calling SNPs in HLA genes by BCFtools call
4.  \-g/--hla-genes, default="all": 26 HLA genes is available in the **hla\_spark**, including A,B,C,DMA,DMB,DOA,DOB,DPA1,DPB1,DQA1,DQA2,DQB1,DQB2,DRA,DRB1,DRB3,DRB4,DRB5,E,F,G,HFE,MICA,MICB,TAP1,TAP2

### Demo for **hla\_spark**&#x20;

Calling HLA SNPs from individual next-generation sequencing (NGS) data

    cd demo
    hla_spark -s sample -f1 fastq/test.r1.fq.gz -f2 fastq/test.r2.fq.gz -n 12 -g A,B,C,DRB1,DQB1

### **hla\_spark** outputs

Compression VCF file containing the HLA SNPs(GRCh38.p14-based)： demo/hla_spark_results/sample/snp_report/**sample.snp.vcf.gz**.&#x20;

The VCF file contains all genomic position with at least one supporting reads of high-quality alignment. &#x20;

