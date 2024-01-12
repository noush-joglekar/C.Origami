import unittest
import corigami.inference.editing as edit
import numpy as np

# Script for testing simple junction prediction code on mskilab server 
#

def testJunction():
	output_path = '/gpfs/commons/home/ieshghi/public_html/contact_junction_background/corigami_testing/'
	celltype = 'imr90' #model claims to be transferable, so let's try it
	
	# let's pick a junction we know exists in H526
	# chr8:60180773 - chr22:30696458
	# left segment is -1 orientation (goes towards increasing genomic coordinates)
	# right segment is +1 orientation 
	
	left_segm_coords = [8,-1,60180773]
	right_segm_coords = [22,+1,30696458]
	
	project_folder = '/gpfs/commons/groups/imielinski_lab/home/ieshghi/Projects/sv_correction_origami/'
	model_path = project_folder + 'corigami_data/model_weights/corigami_base.ckpt'
	seq_path = project_folder + 'corigami_data/data/hg38/dna_sequence'
	ctcf_path = project_folder + 'corigami_data/data/hg38/imr90/genomic_features/ctcf_log2fc.bw'
	atac_path = project_folder + 'corigami_data/data/hg38/imr90/genomic_features/atac.bw'
	
	pred,seq,atac,ctcf = edit.simple_junction(output_path, celltype, left_segm_coords, right_segm_coords, model_path, seq_path, ctcf_path, atac_path)
