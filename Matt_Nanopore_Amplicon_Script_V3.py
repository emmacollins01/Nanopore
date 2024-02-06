#! /usr/bin/env python

'''
=================
Date: 25/07/2023
=================

=================
Prior Processing 
=================

Used Updated Basecalling Script so that it trims nanopore adapters as
filters reads with quality lower than 7 and then  demutltiplexes based on 
Nanopore Barcodes which are also removed. 

= mnt/storageG1/dan/minion/scripts/pipeline_v1.4_demultiplex_emma.sh

Demultiplexed, using demux_nanopore_plates.py script to seperate out
samples based on our barcodes.

Demux_nanopore_plates.py

=================
Assumptions:
=================


1)  Aligner used is BWA MEM as its amplicon data not typical long read data so we dont 
    have to use Minimap2. 


# Temp Delete Me
Aedes_aegypti_lvpagwg_ref.fa = Reference but will have to make better.

python Matt_Nanopore_Amplicon_Script.py --index-file Samples_Index_Subset.csv --ref ~/WGS/Reference_WGS/Aedes_aegypti_lvpagwg_ref.fa --bed Dummy.bed  --ClippingThreshold 50


#################################################
## Command Running to Test Minimap2 (04/09/23) ##
#################################################

 python Matt_Nanopore_Amplicon_Script_V3.py --RunID MINIMAP2Test --index-file ./Sample_Index_Full.csv --ref /mnt/storage8//emma/WGS/Reference_WGS/Aedes_aegypti_lvpagwg_ref.fa --bed Amplicon_Sites.bed --MinReadLength 100 --MaxReadLength 5000 --MinReadQuality 10 --ClippingThreshold 50 --FBMinBaseQual 30 --FBMinCov 50 --gff /mnt/storage8/emma/WGS/Reference_WGS/Aedes_aegypti_lvpagwg.AaegL5.47.gff3 \



'''

import os
import sys
import argparse
import subprocess as sp
import pandas as pd
import csv
import fastq2matrix as fm
from collections import defaultdict
import gzip

# Function to run Bash Commands.
def run_cmd(cmd):
    sys.stderr.write("Running command:\n%s\n\n" % cmd)
    with open("/dev/null","w") as O:
        res = sp.call(cmd,shell=True,stderr=O,stdout=O)
    if res!=0:
        sys.exit("Error running last command, please check!\n")


# Main Function to Run Script
def main(args):

    # Create Run Log
    Log = open(f'{args.RunID}.log','w')
    Log.write(f'\nParameters\n')
    for Param, PV in vars(args).items():
        Log.write(f'{Param}:\t{PV}\n')
    Log.close()


    # Create BWA Index for Reference Genome.
    fm.bwa_index(args.ref)
    fm.create_seq_dict(args.ref)
    fm.faidx(args.ref)

    # Parse Samples IDs from Index File

    # Read in Index CSV File 
    IndexFile = pd.read_csv(args.index_file)
    
    # Ensure only 2 columns of interest present
    IndexFile = IndexFile[['sample','reads']]

    # Drop any rows with null values 
    IndexFile = IndexFile.dropna()

    # Convert 2 columns to dictionary
    RunGuide = dict(zip(IndexFile['sample'], IndexFile['reads']))

    # Create Bam Guide
    BamGuide = open("bam_list.txt","w")

    # Loop over samples. 
    for sample, fastq in RunGuide.items():

        # Update Bam Guide
        BamGuide.write(f'{sample}.bam\n')

        # Check if Bam Exists and Skip if previously made
        if os.path.isfile(f'{sample}.bam')==True:
            print(f'{sample} BAM File Already Exists Skipping!')
            continue

        # Update Arguments
        args.sample = sample
        args.fastq = fastq

        # QC Samples (Write Summary for each Sample) 
        # FiltLong (Link: https://github.com/rrwick/Filtlong)
        #   - Minimum Length 
        #   - Minimum Quality Score. 

        run_cmd(f'filtlong --min_length {args.MinReadLength} --max_length {args.MaxReadLength} --min_mean_q {args.MinReadQuality} {fastq} > {sample}.filtered.fastq')


        # Map Nanopore reads using BWA (LINK:)  
        # Removed Hard and Soft Clipped Ends from read to prevent bad basecalling. (LINK: https://github.com/tseemann/samclip)
        # Sort mapped reads by position using Samtools.
        
        #run_cmd("bwa mem -t 20 -R \"@RG\\tID:%(sample)s\\tSM:%(sample)s\\tPL:nanopore\" %(ref)s %(sample)s.filtered.fastq | samclip --ref %(ref)s --max %(ClippingThreshold)s | samtools sort -o %(sample)s.bam -" % vars(args))
        
        run_cmd("minimap2 -t 20 -a -R \"@RG\\tID:%(sample)s\\tSM:%(sample)s\\tPL:nanopore\" %(ref)s %(sample)s.filtered.fastq | samclip --ref %(ref)s --max %(ClippingThreshold)s | samtools sort -o %(sample)s.bam -" % vars(args))
        
        run_cmd("samtools index %(sample)s.bam" % vars(args))

        
        # Generate Mapping Statistics using Samtools. 
        run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.flagstat.txt" % vars(args))
        
        # Extract Depth Positions and get a summary of the number of positions 
        # in the bed file specific that are covered with at least 1x, 10x, 20x etc
        # (LINK: https://github.com/brentp/mosdepth)
        run_cmd("mosdepth -x -b %(bed)s %(sample)s --thresholds 1,10,20,30  %(sample)s.bam" % vars(args))


        # Calculate the mean coverage for each region specified in the Bed file 
        # (LINK: https://bedtools.readthedocs.io/en/latest/)
        run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_coverage_mean.txt" % vars(args))


    # Close Bam Guide
    BamGuide.close()


    # Call Variants using FreeBayes.(LINK: https://github.com/freebayes/freebayes)
    # Subset VCF based on regions in BED file, noramlise and sort using bcftools.  
    
    if os.path.isfile("MINIMAP2Test.combined.unfiltered.genotyped.vcf.gz")==False:

        run_cmd("freebayes -f %(ref)s -t %(bed)s -L bam_list.txt --haplotype-length -1 --min-coverage %(FBMinCov)s --min-base-quality %(FBMinBaseQual)s --gvcf --gvcf-dont-use-chunk true | bcftools view -T %(bed)s | bcftools norm -f %(ref)s | bcftools sort -Oz -o %(RunID)s.combined.unfiltered.genotyped.vcf.gz" % vars(args))

    print('FILTERING')

    # Filter Variants Based on Depth (10)
    run_cmd("bcftools filter -i 'FMT/DP>10' -S . %(RunID)s.combined.unfiltered.genotyped.vcf.gz | bcftools sort -Oz -o %(RunID)s.combined.DepthFilt.vcf.gz" % vars(args))
    
    # Consequence call SNPs using csq
    run_cmd("bcftools view -v snps %(RunID)s.combined.DepthFilt.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o %(RunID)s.combined.DepthFilt.SNPs.vcf.gz" % vars(args))
    
    # Index Consequence SNPs vcf
    run_cmd("tabix %(RunID)s.combined.DepthFilt.SNPs.vcf.gz" % vars(args))

    # Consequence call Indels. 
    run_cmd("bcftools view -v indels %(RunID)s.combined.DepthFilt.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o %(RunID)s.combined.DepthFilt.INDELs.vcf.gz" % vars(args))
    
    # Index VCF
    run_cmd("tabix %(RunID)s.combined.DepthFilt.INDELs.vcf.gz" % vars(args))
    
    # Convert to Tab-deliminated format. 
    run_cmd("bcftools query {0}.combined.DepthFilt.SNPs.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\n]' > {0}_combined_genotyped_filtered_formatted.snps.txt".format(args.RunID))    
    run_cmd("bcftools query {0}.combined.DepthFilt.SNPs.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\t%TBCSQ\n]' > {0}_combined_genotyped_filtered_formatted.snps.trans.txt".format(args.RunID))
    run_cmd("bcftools query {0}.combined.DepthFilt.INDELs.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\n]' > {0}_combined_genotyped_filtered_formatted.indels.txt".format(args.RunID))
    run_cmd("bcftools query {0}.combined.DepthFilt.INDELs.vcf.gz -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%TGT\t%DP\t%AD\t%TBCSQ\n]' > {0}_combined_genotyped_filtered_formatted.indels.trans.txt".format(args.RunID))



########################################################################################
########################################################################################
########################################################################################
########################################################################################


# Set up the parser
parser = argparse.ArgumentParser(description='Matt Nanopore Amplicon Analysis Script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Essential Arguments
parser.add_argument('--RunID',type=str,help='Run ID used to name log file and VCF',required=True)

parser.add_argument('--index-file',type=str,help='samples.csv with the "sample" column for sample IDs; created by "demux_nanopore_plates.py" if using 96-well plates',required=True)
parser.add_argument('--ref',type=str,help='Reference fasta',required=True)
parser.add_argument('--bed',type=str,help='BED file with genes/amplicon locations',required=True)

# Fastq Filtering Parameters
parser.add_argument('--MinReadLength',type=str,default='100',help='Min Read Length Parameter for Filtlong')
parser.add_argument('--MaxReadLength',type=str,default='5000',help='Max Read Length Parameter for Filtlong')
parser.add_argument('--MinReadQuality',type=str,default='10',help='Min Read Quality Parameter for Filtlong')

# BAM Filtering Parameters used by SamClip 
parser.add_argument('--ClippingThreshold',type=str,default='50',help='Clipping threshold used by samclip')

# FreeBayes Parameters
parser.add_argument('--FBMinBaseQual',default='30',type=int,help='Minimum base quality to use by freebayes')
parser.add_argument('--FBMinCov',default='50',type=int,help='Minimum Coverage to use by freebayes')


# Bcftools consequence calling parameters
parser.add_argument('--gff',type=str,help='GFF')

# Misc Parameters
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.set_defaults(func=main)

# Trigger Script
args = parser.parse_args()
args.func(args)


