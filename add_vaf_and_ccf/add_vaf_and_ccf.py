#!/usr/bin/env python

# This script calculates the variant allele fraction and cancer cell fraction
# of a patients using a MAF file and specific output from TITAN

import argparse
from pybedtools import BedTool
import pandas as pd
import os
import random
import sys
import numpy

def vaf_to_ccf(vaf,pur,prev,minor_cn,major_cn):
    cn = major_cn + minor_cn
    alpha = (cn*prev + 2*(1-prev))*pur + 2*(1-pur)

    #no duplications
    if (minor_cn <= 1 & major_cn <= 1):
        return alpha*vaf/pur

    #One Duplication, no LOH
    elif (minor_cn == 1):
        if (vaf >= (major_cn*prev + 1)*pur/2*alpha):
            return (alpha*vaf/pur)-(prev*(major_cn-1))
        else:
            return alpha*vaf/pur

    #One duplication with LOH
    elif(minor_cn==0):
        if(vaf >= (1+(major_cn-1)*prev)*pur/(2*alpha)):
            return alpha*vaf/pur - (major_cn-1)*prev
        else:
            return alpha*vaf/pur

    #two duplications
    else:
        if(vaf <= (1+(minor_cn-1)*prev)*pur/(2*alpha)):
            return alpha*vaf/pur
        elif(vaf >= (1+(major_cn+minor_cn-1)*prev)*pur/(2*alpha)):
            return alpha*vaf/pur - (major_cn-1)*prev
        else:
            return alpha*vaf/pur - (minor_cn-1)*prev

if __name__ == "__main__":
    desc = ""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--replace", "-r", action="store_true", default=False)
    parser.add_argument("--distance", "-d", type=int, nargs=1, default=3200000000)
    parser.add_argument("maf")
    parser.add_argument("params")
    parser.add_argument("segs")
    parser.add_argument("outfile")
    args = parser.parse_args()


    ### GRAB THE PURITY AND VALIDITY INDEX FROM PARAMS.TXT ############
    param = pd.io.parsers.read_table(
        args.params,
        comment='#',
        header=None,
        index_col=0,
        sep=':\t'
        )

    purity = 1 - float(param.loc[['Normal contamination estimate'],1][0])
    validity = param.loc[['S_Dbw validity index'],1][0]


    ### PARSE MAF FILE AND APPEND NEW COLUMNS #########################
    maf = pd.io.parsers.read_table(
        args.maf,
        header=0,
        comment='#',
        sep='\t'
        )

    # Required Columns (keys) in the MAF file
    count_keys = ["n_ref_count", "n_alt_count", "t_ref_count", "t_alt_count"]
    for key in count_keys:
        if not key in maf.dtypes:
            sys.exit(''.join(["The input MAF must contain the following column: ", key]))

    # Grab the dimensions of the data frame
    nrow_maf = maf.shape[0]
    ncol_maf = maf.shape[1]

    # New Columns (keys) to add to the MAF file if not already present
    new_keys = ['purity', 'prevelence', 'minor_cn', 'major_cn', 'vaf', 'ccf']
    for key in new_keys:
        if not key in maf.dtypes or args.replace:
            maf[key] = pd.Series([0.0] * nrow_maf, index=maf.index)


    ### PARSE SEG FILE ################################################
    seg = pd.io.parsers.read_table(
        args.segs,
        header=0,
        comment='#',        
        sep='\t'
        )

    nrow_seg = seg.shape[0]
    ncol_sef = seg.shape[1]


    ### CREATE SUB BED FILES FROM MAF AND SEG #########################
    maf_sub_dataframe = maf[['Chromosome', 'Start_Position', 'End_Position']]
    maf_sub_dataframe['Maf_index'] = range(0, nrow_maf)
    seg_sub_dataframe = seg[['Chromosome', 'Start_Position(bp)','End_Position(bp)']]
    seg_sub_dataframe['Seg_index'] = range(0, nrow_seg)

    maf_bed = './maf.bed'
    seg_bed = './seg.bed'

    while os.path.isfile(maf_bed):
        maf_bed = '.'.join(['./maf', str(random.randint(100000,999999)), 'bed'])
    while os.path.isfile(seg_bed):
        seg_bed = '.'.join(['./seg', str(random.randint(100000,999999)), 'bed'])

    maf_sub_dataframe.to_csv(maf_bed, sep='\t', header=0, index=False)
    seg_sub_dataframe.to_csv(seg_bed, sep='\t', header=0, index=False)

    files_to_remove = [maf_bed, seg_bed]
    maf_bed = BedTool(maf_bed)
    maf_bed = maf_bed.sort()
    seg_bed = BedTool(seg_bed)
    seg_bed = seg_bed.sort()

    closest_bed = maf_bed.closest(seg_bed, d=True)

    closest_dataframe = closest_bed.to_dataframe()
    closest_dataframe.columns = ['Maf_Chrom', 'Maf_Start', 'Maf_End', 'Maf_Index', 'Seg_Chrom', 'Seg_Start', 'Seg_End', 'Seg_Index', 'Distance']

    nrow_closest = closest_dataframe.shape[0]
    ncol_closest = closest_dataframe.shape[1]


    ### CALCULATE VAF & CCF ###########################################
    # Pulls together all information needed for calculation and adds
    # data to the appropriate column in maf by using indexes which
    # were stored in the maf_bed and seg_bed tmp bed data types
    # Allows for O(n) runtime, without having to serach and find
    # the specific index in the MAF where the chr, start and end match
    for row in range(0, nrow_closest):

        if closest_dataframe['Distance'] > args.distance:
            continue
        
        # Maf and Seg Indexes
        maf_index = closest_dataframe['Maf_Index'][row]
        seg_index = closest_dataframe['Seg_Index'][row]
        
        # Some chromosomes which appear in the MAF are not present in 
        # The segs file, these are represented with a '.', so if
        # This is such a region, continue to the next line
        if seg_index == '.' or maf_index == '.':
            continue
        
        # Otherwise maf_index and seg_index is an integer and should
        # bed converted to their appropriate type
        else:
            maf_index = int(maf_index)
            seg_index = int(seg_index)

        # Pull all the variables from the appropriate files
        t_alt_count = float(numpy.int(maf['t_alt_count'][maf_index]))
        t_ref_count = float(numpy.int(maf['t_ref_count'][maf_index]))
        prev = float(seg['Clonal_Frequency'][seg_index])
        minor_cn = seg['MinorCN'][seg_index]
        major_cn = seg['MajorCN'][seg_index]
        
        # Calculate VAF
        vaf = float(t_alt_count / (t_alt_count + t_ref_count))

        # Calcualte CCF
        ccf = float(vaf_to_ccf(vaf, purity, prev, minor_cn, major_cn))
        
        # Store in maf data type
        maf['vaf'][maf_index] = float(vaf)
        maf['purity'][maf_index] = float(purity)
        maf['prevelence'][maf_index] = float(prev)
        maf['minor_cn'][maf_index] = float(minor_cn)
        maf['major_cn'][maf_index] = float(major_cn)
        maf['ccf'][maf_index] = float(ccf)


    ### WRITE OUTFILE #################################################
    maf.to_csv(
        args.outfile,
        sep = '\t',
        header = True,
        index = False
        )


    ### CLEAN TEMP FILES ##############################################
    os.remove(files_to_remove[0])
    os.remove(files_to_remove[1])

