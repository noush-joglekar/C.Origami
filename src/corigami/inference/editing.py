import os
import numpy as np
import pandas as pd
import sys
import torch

import corigami.inference.utils.inference_utils as infer
from corigami.inference.utils import plot_utils 

import argparse

def main():
    parser = argparse.ArgumentParser(description='C.Origami Editing Module.')
    
    # Output location
    parser.add_argument('--out', dest='output_path', 
                        default='outputs',
                        help='output path for storing results (default: %(default)s)')

    # Location related params
    parser.add_argument('--celltype', dest='celltype', 
                        help='Sample cell type for prediction, used for output separation', required=True)
    parser.add_argument('--chr', dest='chr_name', 
                        help='Chromosome for prediction', required=True)
    parser.add_argument('--start', dest='start', type=int,
                        help='Starting point for prediction (width is 2097152 bp which is the input window size)', required=True)
    parser.add_argument('--model', dest='model_path', 
                        help='Path to the model checkpoint', required=True)
    parser.add_argument('--seq', dest='seq_path', 
                        help='Path to the folder where the sequence .fa.gz files are stored', required=True)
    parser.add_argument('--ctcf', dest='ctcf_path', 
                        help='Path to the folder where the CTCF ChIP-seq .bw files are stored', required=True)
    parser.add_argument('--atac', dest='atac_path', 
                        help='Path to the folder where the ATAC-seq .bw files are stored', required=True)

    # Deletion related params.
    # Noush - making them generalizable for SVs
    parser.add_argument('--svBlock-start', dest='sv_start', type=int,
                        help='Starting point for deletion or insertion block on diff chrom.', required=True)
    parser.add_argument('--svBlock-width', dest='sv_width', type=int,
                        help='Width for deletion or insertion', required=True)
    parser.add_argument('--padding', dest='end_padding_type', 
                        default='zero',
                        help='Padding type, either zero or follow. Using zero: the missing region at the end will be padded with zero for ctcf and atac seq, while sequence will be padded with N (unknown necleotide). Using follow: the end will be padded with features in the following region (default: %(default)s)')
    parser.add_argument('--hide-line', dest='hide_deletion_line', 
                        action = 'store_true',
                        help='Remove the line showing deletion site (default: %(default)s)')

    parser.add_argument('--sv-type', dest='sv_type', type=str, default="deletion")
    # Insertion related params
    parser.add_argument('--insertion-source-start',dest='ins_site', type=int) ### SWITCH THIS WITH DELETION_START
    parser.add_argument('--insertion-source-chrom', dest='ins_chrom')

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if args.sv_type.lower() == "deletion":
        single_deletion(args.output_path, args.celltype, args.chr_name, args.start, 
                        args.sv_start, args.sv_width, 
                        args.model_path,
                        args.seq_path, args.ctcf_path, args.atac_path, 
                        show_deletion_line = not args.hide_deletion_line,
                        end_padding_type = args.end_padding_type)
    if args.sv_type.lower() == "insertion":
        print("Insertion selected ... ")
        single_insertion(args.output_path, args.celltype, args.chr_name, args.start, 
                        args.sv_start, args.sv_width, 
                        args.model_path,
                        args.seq_path, args.ctcf_path, args.atac_path,
                        args.ins_site, args.ins_chrom,
                        show_insertion_lines = not args.hide_deletion_line,
                        end_padding_type = args.end_padding_type)

def single_deletion(output_path, celltype, chr_name, start, deletion_start, deletion_width, model_path, seq_path, ctcf_path, atac_path, show_deletion_line = True, end_padding_type = 'zero'):

    # Define window which accomodates deletion
    window = 2097152 + deletion_width
    seq_region, ctcf_region, atac_region = infer.load_region(chr_name, 
            start, seq_path, ctcf_path, atac_path, window = window)
    # Delete inputs
    seq_region, ctcf_region, atac_region = deletion_with_padding(start, 
            deletion_start, deletion_width, seq_region, ctcf_region, 
            atac_region, end_padding_type)
    # Prediction
    pred = infer.prediction(seq_region, ctcf_region, atac_region, model_path)
    # Initialize plotting class
    plot = plot_utils.MatrixPlotDeletion(output_path, pred, 'deletion', 
            celltype, chr_name, start, deletion_start, deletion_width, 
            padding_type = end_padding_type,
            show_deletion_line = show_deletion_line)
    plot.plot()

def deletion_with_padding(start, deletion_start, deletion_width, seq_region, ctcf_region, atac_region, end_padding_type):
    ''' Delete all signals at a specfied location with corresponding padding at the end '''
    # end_padding_type takes values of either 'zero' or 'follow'
    if end_padding_type == 'zero':
        seq_region, ctcf_region, atac_region = zero_region(seq_region, 
                ctcf_region, atac_region)
    elif end_padding_type == 'follow': pass
    else: raise Exception('unknown padding')
    # Deletion
    seq_region, ctcf_region, atac_region = delete(deletion_start - start, 
            deletion_start - start + deletion_width, 
            seq_region, ctcf_region, atac_region)
    return seq_region, ctcf_region, atac_region

def zero_region(seq, ctcf, atac, window = 2097152):
    ''' Replace signal with zero. N for sequence and 0 for CTCF and ATAC '''
    seq[window:] = [0, 0, 0, 0, 1]
    ctcf[window:] = 0
    atac[window:] = 0
    return seq, ctcf, atac

def delete(start, end, seq, ctcf, atac, window = 2097152):
    seq = np.delete(seq, np.s_[start:end], axis = 0)
    ctcf = np.delete(ctcf, np.s_[start:end])
    atac = np.delete(atac, np.s_[start:end])
    return seq[:window], ctcf[:window], atac[:window]


def single_insertion(output_path, celltype, chr_name, start, insertion_start, insertion_width, model_path, seq_path, ctcf_path, atac_path, 
ins_site, ins_chrom, show_insertion_lines = True, end_padding_type = 'zero'):
    """Replace deletion with insertion. Just add an extra parameter detailing which function to call. Should be equivalent otherwise"""
    window = 2097152 - insertion_width ## carrying over. i think it has to do with bits.
    hWindow = window // 2
    left_edge = insertion_start - (insertion_width // 2 ) - hWindow ## not sure if this works. we'll see
    right_edge = insertion_start + (insertion_width // 2 ) + hWindow
    ins_seq, ins_ctcf, ins_atac = infer.load_region(ins_chrom, ins_site, seq_path, ctcf_path, atac_path, window = insertion_width)
    left_seq, left_ctcf, left_atac = infer.load_region(chr_name, left_edge, seq_path, ctcf_path, atac_path, window=hWindow)
    right_seq, right_ctcf, right_atac = infer.load_region(chr_name, right_edge, seq_path, ctcf_path, atac_path, window=hWindow)

    ## join the inputs:
    seq_region = np.concatenate((left_seq,ins_seq,right_seq))
    atac_region = np.concatenate((left_atac,ins_atac,right_atac))
    ctcf_region = np.concatenate((left_ctcf,ins_ctcf,right_ctcf))

    # Prediction
    pred = infer.prediction(seq_region, ctcf_region, atac_region, model_path)
    print("Insertion stuff done ... plotting now")
    # # Initialize plotting class
    plot = plot_utils.MatrixPlotInsertion(output_path, pred, 'insertion', 
            celltype, chr_name, start, insertion_start, insertion_width, ins_chrom,
            padding_type = end_padding_type,
            show_insertion_lines = show_insertion_lines)
    
    plot.plot()
    print("Plotted ??")

# A simple junction between two chromosomes
# Added by: Iraj , m-ski lab
#

def simple_junction(output_path, celltype, left_segm_coords, right_segm_coords, model_path, seq_path, ctcf_path, atac_path, show_deletion_line = True, bp_loc_in_window = None, ifplot = True):
# new parameters are:
#	left_segm_coords = [chromosome number, orientation, breakpoint] 
#	right_segm_coords = [chromosome number, orientation, breakpoint] 
# 	bp_loc_in_window: location of the breakpoint in the simulation window. If None, the window is centered at the breakpoint. Otherwise, a floating-point number between 0 and 1 gives the fractional location along the window
	if bp_loc_in_window is None:
		bp_loc_in_window = 0.5

	assert isinstance(bp_loc_in_window,float) and bp_loc_in_window >=0 and bp_loc_in_window <=1
# 
# Orientations are either +1 or -1
# For the "left" segment: if the orientation is +1, the segment goes from coordinate 1 to breakpoint. Otherwise, the segment goes from breakpoint to the "end" of the chromosome
# For the "right" segment: if the orientation is +1, the segment goes from breakpoint to the "end" of the chromosome. Otherwise, the segment goes from coordinate 1 to the breakpoint
#
# 
# The other parameters function identically as in the function single_deletion() above

	left_chr_name = 'chr'+ str(left_segm_coords[0])
	right_chr_name = 'chr'+ str(right_segm_coords[0])

    # Define window which accomodates SV
	window = 2097152 #this is hard-coded into the model...? might need to add a layer if using a different window size
	left_length = int(bp_loc_in_window*window)
	right_length = window - left_length #guarantees that the resulting snippet of sequence is the right length

	# get the sequences for the left side of the breakpoint, flip if necessary
	left_breakpoint = left_segm_coords[2]
	right_breakpoint = right_segm_coords[2]
	assert (left_segm_coords[1] in [1,-1]) and (right_segm_coords[1] in [1,-1])
	if left_segm_coords[1] == 1: #if orientation = +1, left segment spans [breakpoint - left_length , breakpoint]
		seq_left, ctcf_left, atac_left = infer.load_region(left_chr_name, left_breakpoint - left_length, seq_path, ctcf_path, atac_path, window = left_length)
	else: #if orientation is -1, left segment goes [breakpoint, breakpoint + left_length] and is then inverted
		seq_left, ctcf_left, atac_left = infer.load_region(left_chr_name, left_breakpoint, seq_path, ctcf_path, atac_path, window = left_length)
		#now we must invert the sequences
		seq_left = np.flip(seq_left,0).copy()
		atac_left = np.flip(atac_left,0).copy()
		ctcf_left = np.flip(ctcf_left,0).copy()

	# same but opposite on the right side
	if right_segm_coords[1] == 1:#if orientation = +1, right segment spans [breakpoint, breakpoint + right_length]
		seq_right, ctcf_right, atac_right = infer.load_region(right_chr_name, right_breakpoint, seq_path, ctcf_path, atac_path, window = right_length)
	else: #if orientation = -1, right segment spans [breakpoint-right_length, breakpoint] and is then inverted
		seq_right, ctcf_right, atac_right = infer.load_region(right_chr_name, right_breakpoint - right_length, seq_path, ctcf_path, atac_path, window = right_length)
		#now we invert the sequences
		seq_right = np.flip(seq_right,0).copy()
		atac_right = np.flip(atac_right,0).copy()
		ctcf_right = np.flip(ctcf_right,0).copy()

	#splice the inputs
	seq_region = np.concatenate((seq_left,seq_right))
	atac_region = np.concatenate((atac_left,atac_right))
	ctcf_region = np.concatenate((ctcf_left,ctcf_right))
		
    	#prediction
	pred = infer.prediction(seq_region, ctcf_region, atac_region, model_path)

#	plot = plot_utils.MatrixPlot(output_path, pred, 'simple_junction', celltype, 
#                                 left_chr_name + '-' + right_chr_name, left_breakpoint-left_length)
	#we made our own plotting function. does it work?

	if ifplot:
		plot = plot_utils.MatrixPlotJunction(output_path, pred, 'simple_junction', celltype,left_segm_coords,right_segm_coords,bp_loc_in_window)
		plot.plot()

	return pred

if __name__ == '__main__':
    main()

