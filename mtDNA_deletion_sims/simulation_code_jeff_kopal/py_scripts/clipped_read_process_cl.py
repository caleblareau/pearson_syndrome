import re
import os
import pysam
from collections import Counter
import sys
from optparse import OptionParser

# Parse out data
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files and split based on some attribute"
opts = OptionParser(usage=usage)
opts.add_option("-i", "--input", help="Filename of new .bam file to be generated")
opts.add_option("-b", "--hq-barcodes-file", help="Filepath of the high quality barcodes", 
                default=None)
opts.add_option("-t", "--sam-tag", help="two letter sam tag needed for file splitting")

options, arguments = opts.parse_args()

inbam = options.input
hq_barcodes_file = options.hq_barcodes_file
tagofinterest = options.sam_tag

bam_in = pysam.AlignmentFile(inbam, "rb")

    
def process_cigar_for_clip_position(cigar, tuple):

    pos = 0

    # Case 1/2: start of read
    if(cigar[1] == "H" or cigar[2] == "H" or
       cigar[1] == "S" or cigar[2] == "S"):
        pos = tuple[0][1] + 1 #offset 1-base and differentiate 0 pos/noclip

    # Case 2: end of read
    if(cigar[-1] == "H" or cigar[-1] == "S"):
        pos = tuple[-1][1] + 1 #offset 1-base and differentiate 0 pos/noclip

    return(pos)
        
def left_clip(cigar):
    if(cigar[1] == "H" or cigar[2] == "H" or
       cigar[1] == "S" or cigar[2] == "S"):
        return(1)
    else:
        return(0)

def right_clip(cigar):
    if(cigar[-1] == "S" or cigar[-1] == "H"):
        return(1)
    else:
        return(0)

def getCBvalue(input_tags):
    for tg in input_tags:
        if("CB" == tg[0]):
            return(tg[1])
    return("NA")
    
def getTag(intags, tag):
    '''
    Checks for specific tag from CLI
    '''
    for tg in intags:
        if(tag == tg[0]):
            return(tg[1])
    return("NA")




# Open all the output files and spit out the filtered data
# Based on where the matching value originates

clip_pos_count = Counter()
outfile = re.sub(".bam$", ".stats.tsv", inbam)
outfile_handle = open(outfile, 'w')

if(hq_barcodes_file is None):
    for read in bam_in.fetch('chrM'):
        seq = read.seq
        reverse = read.is_reverse
        #quality = read.query_qualities
        #align_qual_read = read.mapping_quality
        cigar_string = read.cigarstring
        positions = read.get_reference_positions()
        tuple = read.get_aligned_pairs(True)
        if(positions and cigar_string):
            start = str(positions[0] + 1) #offset to 1 base
            end = str(positions[-1] + 1) #offset to 1 base
            #rev.append(reverse)
            rc = str(right_clip(cigar_string))
            lc = str(left_clip(cigar_string))
            cell_id = str(getCBvalue(read.tags))
            clip_pos = str(process_cigar_for_clip_position(cigar_string, tuple))
            clip_pos_count[clip_pos] += 1
            read_name = str(read.query_name)
            list_of_outs = [start, end, lc, rc, cell_id, clip_pos, read_name]
            outfile_handle.write("\t".join(list_of_outs) + "\n")
else:
    # Parse a dictionary file
    hq_barcodes = set([line.rstrip('\n') for line in open(hq_barcodes_file)])
    for read in bam_in.fetch('chrM'):
        tag = getTag(read.tags, tagofinterest)
        if(tag in hq_barcodes):
            seq = read.seq
            reverse = read.is_reverse
            #quality = read.query_qualities
            #align_qual_read = read.mapping_quality
            cigar_string = read.cigarstring
            positions = read.get_reference_positions()
            tuple = read.get_aligned_pairs(True)
            if(positions and cigar_string):
                start = str(positions[0] + 1) #offset to 1 base
                end = str(positions[-1] + 1) #offset to 1 base
                #rev.append(reverse)
                rc = str(right_clip(cigar_string))
                lc = str(left_clip(cigar_string))
                cell_id = str(getCBvalue(read.tags))
                clip_pos = str(process_cigar_for_clip_position(cigar_string, tuple))
                clip_pos_count[clip_pos] += 1
                read_name = str(read.query_name)
                list_of_outs = [start, end, lc, rc, cell_id, clip_pos, read_name]
                outfile_handle.write("\t".join(list_of_outs) + "\n")
outfile_handle.close()

#write out clip_pos counter as table
outfile_clip = re.sub(".bam$", ".clip.tsv", inbam)
with open(outfile_clip, mode='w') as fp:
    fp.write('clip_pos\tcount\n')
    for tag, count in clip_pos_count.items():
        fp.write('{}\t{}\n'.format(tag, count)) 


