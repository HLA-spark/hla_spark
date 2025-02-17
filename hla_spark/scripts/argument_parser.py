import argparse
import os

usage = '''
HLA-spark is a tool for calling coding SNPs in HLA genes from NGS-data. 
    Results (Human genomic reference: GRCh38.p14): 
    1) Compression VCF file containing BCFtools called SNPs in HLA genes-target sites : analysis_path/$SAMPLE/snp_report/$SAMPLE.snp.vcf.gz.
    2) Compression VCF file containing genotype information that can be piped to BCFtools to call HLA SNPs:  analysis_path/$SAMPLE/snp_report/$SAMPLE.genotype.vcf.gz.
Example: 
    python  %s  -h
''' % (__file__[__file__.rfind(os.sep) + 1:])

def argument_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=usage)

    #Required arguments
    required=parser.add_argument_group('Required arguments')
    required.add_argument('-s','--sample_name',help='input the sample name, used for naming the analysis directory',
                          type=str,default=None,required=True)
    required.add_argument('-f1', '--fastq1', help='input the file path of r1.fq(.gz) data',
                         type=str, default=None,required=True)
    required.add_argument('-f2', '--fastq2', help='input the file path of r2.fq(.gz) data',
                          type=str, default=None, required=True)

    #Options
    parser.add_argument('-o', '--output-dir', type=str, default="hla_spark_results/",
                        help="Directory path for output. You can put the analysis results of multiple samples under a same directory. A subdirectory named the input sample name will be auto-created under the output-dir.\ndefault='hla_spark_results/'")
    parser.add_argument('-n', '--n-core', type=int, default=4,
                        help="Number of threads to use in the analysis\ndefault = 4")
    parser.add_argument('-l', '--step-list', type=str, default="1,2,3,4,5,6,7,8", help='''Run Steps:
            1 : HLA_reads_extraction_module
            2 : Best_matched_alleles_module
            3 : Align_module
            4 : Rmdup_module
            5 : Mpileup_module
            6 : VCFConversion_module
            7 : Genotype_report_module
            8 : SNP_calling_module
            default = 1,2,3,4,5,6,7,8
        ''')
    parser.add_argument('-g','--hla-genes', type=str, default='all',help='''26 avalible hla genes can be selected. 
            Default = all: analysis on all of genes. 
            Analysis on one (-g A) or multiple genes (-g A,DRB1,DQA1) is permitted. 
            List of 26 genes: A,B,C,DMA,DMB,DOA,DOB,DPA1,DPB1,DQA1,DQA2,DQB1,DQB2,DRA,DRB1,DRB3,DRB4,DRB5,E,F,G,HFE,MICA,MICB,TAP1,TAP2
    ''')
    arguments = parser.parse_args()
    return arguments.__dict__