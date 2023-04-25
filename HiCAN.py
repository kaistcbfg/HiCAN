
import numpy as np
import pandas as pd

from itertools import product

import argparse
import pickle
import gzip
import sys
import os

import warnings
warnings.filterwarnings("ignore")

total_step = 5 

def str2bool(v):
	if isinstance(v, bool): return v
	if v.lower() in ['yes', 'true', 't', 'y', '1']: return True
	elif v.lower() in ['no', 'false', 'f', 'n', '0']: return False
	else: raise argparse.ArgumentTypeError('Boolean value expected.')
#

def gen_genomeBin_list(chrlist_filepath, bindir, genome, resolution):

	f = open(chrlist_filepath)
	chrlist = [i.rstrip() for i in f.readlines()]
	f.close()

	all_bins = []
	for chrname in chrlist:
		binfile = '{}/{}.{}.{}.bin'.format(bindir, genome, chrname, resolution)
		f = open(binfile)
		for line in f:
			line = line.rstrip()
			binid = ".".join(line.split())
			all_bins.append(binid)
		#
		f.close()
	#
	print("chrlist: {}".format(chrlist))
	print("Total {} bins loaded | resolution: {}".format(len(all_bins),resolution))
	print("bin files from {} loaded. Genome version: {}".format(bindir, genome))

	return all_bins
#

def parse_text_imported_Bins(all_bins, resolution):

	chrlist = [i.split('.')[0] for i in all_bins]
	chrlist = list(set(chrlist))
	chrlist.sort()

	print("chrlist: {}".format(chrlist))
	print("Total {} bins loaded | resolution: {}".format(len(all_bins),resolution))	
	print("row/cols from input file imported. --genome-bindir, --genome-version, and --chrname-list options ignored")
	return
#

def parse_blacklist_file(blacklist_file_path, resolution):

	blk_binidlist = []

	f = open(blacklist_file_path)
	for line in f:
		line = line.rstrip()
		[chrname, start, end] = line.split()

		start = int(start)
		bin1 = start//resolution
		binid1 = ".".join([chrname, str(bin1 * resolution), str((bin1 + 1)* resolution)])

		end = int(end)
		bin2 = end//resolution
		binid2 = ".".join([chrname, str(bin2 * resolution), str((bin2 + 1)* resolution)])

		blk_binidlist.append(binid1)
		blk_binidlist.append(binid2)
	#
	f.close()
	blk_binidlist = list(set(blk_binidlist))
	blk_binidlist.sort()
	print("Total {} blacklist bins in resolution of {}".format(len(blk_binidlist), resolution))

	return blk_binidlist
#

def get_topK_bins(Wval_list, K):
	
	maxW_bins_list = []

	for i in range(3):
		Wval_list = sorted(Wval_list, key=lambda tup:tup[i+1], reverse=True)
		Wval_list_top = Wval_list[:K]
		maxW_bins = [i for i,j,k,l in Wval_list_top]
		maxW_bins_list.append(maxW_bins)
	#

	return maxW_bins_list
#

def make_geneBin_dict(annot_file, resolution):

	gene_binDict = {}

	f = open(annot_file)
	for line in f:
		line = line.rstrip()
		linedata = line.split()

		chrname = linedata[0]
		start = int(linedata[3])
		end = int(linedata[4])
		strand = linedata[6]

		geneid = linedata[9][1:-1]
		genename = linedata[17][1:-1]
		geneinfo = (geneid + "_" + genename)

		if strand == '+': tss = start
		else: tss = end
		tss_binidx = tss//resolution

		gene_bin = ".".join([chrname, str(tss_binidx * resolution), str((tss_binidx+1) * resolution)])

		if gene_bin in gene_binDict.keys():
			if not gene_bin in gene_binDict[gene_bin]: gene_binDict[gene_bin].append(geneinfo)
		else: gene_binDict[gene_bin] = [geneinfo]

	#
	f.close()

	return gene_binDict
#

def get_genelist(binlist, gene_binDict):

	outlist = []

	for nmf_bin in binlist:
		genes_in_bin = gene_binDict.get(nmf_bin)
		if genes_in_bin != None: outlist += genes_in_bin
	#
	return outlist
#

def write_basis_file(outfilename, Wval_list, genecountlist):

	o = open(outfilename,'wt')
	o.write('#{}\t{}\t{}\n'.format(genecountlist[0], genecountlist[1], genecountlist[2]))
	for binW in Wval_list:
		outstr = '\t'.join([str(i) for i in binW])
		outstr += '\n'
		o.write(outstr)
	#
	o.close()

	return
#

def write_gene_file(outfilename, geneinfolist):

	o = open(outfilename,'wt')
	for i in geneinfolist:
		outstr = '\t'.join(i.split('_'))
		outstr += '\n'
		o.write(outstr)
	#
	o.close()

	return
#


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='HiCAN: inter-chromosomal hub interaction caller')
	parser.add_argument('--header', type=str, help='save file header', required=True)
	parser.add_argument('--input-file', type=str, help='input *.gz file (required)', required=True)
	parser.add_argument('--input-format', type=str, default='pickle', help='default pickle, format: text or pickle (numpy array pickle)')
	parser.add_argument('--gzip-flag', type=str2bool, default=True, help='default True, *.gz input')
	parser.add_argument('--resolution', type=int, default=500000, help='bin resolution (default 500kb')
	parser.add_argument('--nmf-solver', type=str, default='R', help='default R, format: R or Python')
	parser.add_argument('--randstate-py', type=int, default=None, help='random state value, default None')
	parser.add_argument('--output-dir', type=str, default='./output', help='default ./output')

	parser.add_argument('--genome-version', type=str, default='hg19', help='default hg19, [hg19, hg38, mm10] available')
	parser.add_argument('--genome-bindir', type=str, default='./genome_Bin', help='default ./genome_Bin. Look [genome].[chrname].[resolution].bin file in dir')
	parser.add_argument('--chrname-list', type=str, default='./CHRLIST.txt', help='default ./CHRLIST.txt List of chromosomes')

	parser.add_argument('--coverage-cutoff', type=float, default=-1, help='default 1/3 of mean coverage (-1), should be > 0, row/col > thresh will remain')
	parser.add_argument('--blacklist-file', type=str, default=False, help='default False, remove corresponding rows/cols from input matrix')
	parser.add_argument('--ICE-flag', type=str2bool, default=False, help='default False, apply default ICE norm to matrix')

	parser.add_argument('--annot-file', type=str, help='gene annotation file', required=True)
	parser.add_argument('--topK-basis', type=int, default=500, help='default 500, get top K n-th basis to calc gene density')	

	parser.add_argument('--visualize-flag', type=str2bool, default=True, help='Save contact maps to PDFs')
	parser.add_argument('--visualize-vmax', type=float, default=50, help='default 50, visualize vmax value')
	parser.add_argument('--visualize-output', type=str, default='./PDFs', help='PDF file save path')

	args = parser.parse_args()

	if args.visualize_flag:
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
	#
	if args.nmf_solver == 'R':
		import rpy2.robjects as ro
		from rpy2.robjects import r
		from rpy2.robjects import pandas2ri
		from rpy2.robjects.conversion import localconverter
		from rpy2.robjects.packages import importr
		NMFr = importr('NMF')	
	elif args.nmf_solver == 'Python': from sklearn.decomposition import NMF
	
	if args.input_format == 'pickle':
		if args.gzip_flag:
			with gzip.open(args.input_file, 'rb') as f: rawmat = pickle.load(f)
		else:
			with open(args.input_file, 'rb') as f: rawmat = pickle.load(f)
		#
		all_bins = gen_genomeBin_list(args.chrname_list, args.genome_bindir, args.genome_version, args.resolution)
	elif args.input_format == 'text':
		if args.gzip_flag: df = pd.read_csv(args.input_file, sep='\t', header=0, index_col=0, compression='gzip')	
		else: df = pd.read_csv(args.input_file, sep='\t', header=0, index_col=0)
		rawmat = df.to_numpy()
		rawmat = np.nan_to_num(rawmat)	
		all_bins = df.columns
		parse_text_imported_Bins(all_bins, args.resolution)
	else: sys.exit('Error: invalid --input-format value') 
	
	print("{} processing started....".format(args.header))
	print("[1/{}] Data imported".format(total_step))
	if args.visualize_flag:
		fig = plt.figure(1)
		ax = fig.add_subplot(111)
		cax = ax.matshow(rawmat, cmap=plt.cm.RdYlBu_r, vmin=0, vmax=args.visualize_vmax)
		fig.colorbar(cax)
		plt.savefig('{}/{}_rawMat.pdf'.format(args.visualize_output, args.header), dpi=1000)
	#

	df = pd.DataFrame(rawmat, columns = all_bins, index = all_bins)

	if args.blacklist_file != False: blacklist_bins = parse_blacklist_file(args.blacklist_file, args.resolution) 
	else: blacklist_bins = []
	df_blkfilt = df.loc[~df.index.isin(blacklist_bins) , ~df.columns.isin(blacklist_bins)]

	if float(args.coverage_cutoff) < 0: coverage_thresh = np.mean(np.sum(df_blkfilt.to_numpy(),axis=0))/3
	else: coverage_thresh = args.coverage_cutoff
	df_cov_blkfilt = df_blkfilt.loc[(df_blkfilt.sum(axis=1) > coverage_thresh), (df_blkfilt.sum(axis=0) > coverage_thresh)]

	if args.ICE_flag:
		from iced import normalization
		icemat = normalization.ICE_normalization(df_cov_blkfilt.to_numpy())
		icemat = np.nan_to_num(icemat)
		df_cov_blkfilt = pd.DataFrame(icemat, columns=df_cov_blkfilt.columns, index=df_cov_blkfilt.index)
		print("ICE normalization applied")
	else: print("ICE normalization OFF")

	cismask = [i.split('.')[0]!=j.split('.')[0] for i,j in product(df_cov_blkfilt.columns,df_cov_blkfilt.columns)]
	cismask = np.array(cismask).reshape((len(df_cov_blkfilt.columns),len(df_cov_blkfilt.columns)))
	df_cis_cov_blkfilt = df_cov_blkfilt.where(cismask)
	df_cis_cov_blkfilt = df_cis_cov_blkfilt.fillna(0)
	df_cis_cov_blkfilt = df_cis_cov_blkfilt.loc[(df_cis_cov_blkfilt.sum(axis=1) > 0), (df_cis_cov_blkfilt.sum(axis=0) > 0)]

	print("[2/{}] Input preprocessed. {} rows/cols remain after filtered".format(total_step,len(df_cis_cov_blkfilt.columns)))
	if args.visualize_flag:
		fig = plt.figure(2)
		ax = fig.add_subplot(111)
		cax = ax.matshow(df_cis_cov_blkfilt.to_numpy(), cmap=plt.cm.RdYlBu_r, vmin=0, vmax=args.visualize_vmax)
		fig.colorbar(cax)
		plt.savefig('{}/{}_filteredMat.pdf'.format(args.visualize_output, args.header), dpi=1000)
	#

	if args.nmf_solver == 'R':
		with localconverter(ro.default_converter + pandas2ri.converter): df_r = ro.conversion.py2rpy(df_cis_cov_blkfilt)
		r.assign('df_r',df_r)
		r('nmf_model <- nmf(df_r, 3, seed = "nndsvd")')
		r('nmf_model_fit <- fit(nmf_model)')
		r('nmf_model_fit_w <- nmf_model_fit@W')
		W_r = r('nmf_model_fit_w')
		W = np.asarray(W_r)
	elif args.nmf_solver == 'Python':
		nmf = NMF(init='nndsvd',n_components=3,random_state=args.randstate_py)
		W = nmf.fit_transform(df_cis_cov_blkfilt.to_numpy())
	else:
		print("Invalid solver type: choose 'R' or 'Python'")
		sys.exit()
	#	
	print("[3/{}] {} NMF decomposition applied".format(total_step, args.nmf_solver))

	max_basis_arr = np.argmax(W,axis=1)
	all_basis_index = []
	basis_len_list = []
	for i in range(3):
		basis_val_list = list(zip(df_cis_cov_blkfilt.columns[max_basis_arr == i] ,W[:,i][max_basis_arr == i]))
		basis_val_list = sorted(basis_val_list, key = lambda tup: tup[1], reverse=True)
		basis_len_list.append(len(basis_val_list))
		for i,j in basis_val_list: all_basis_index.append(i)
	#
	print("[4/{}] {} bins in each NMF basis".format(total_step, basis_len_list))

	df_reorder = df_cis_cov_blkfilt.reindex(index=all_basis_index,columns=all_basis_index)
	if args.visualize_flag:
		fig = plt.figure(3)
		ax = fig.add_subplot(111)
		cax = ax.matshow(df_reorder.to_numpy(), cmap=plt.cm.RdYlBu_r, vmin=0, vmax=args.visualize_vmax)
		fig.colorbar(cax)
		plt.savefig('{}/{}_reorderedMat.pdf'.format(args.visualize_output, args.header), dpi=1000)
	#

	Wval_list = list(zip(df_cis_cov_blkfilt.columns,W[:,0],W[:,1],W[:,2]))
	gene_binDict = make_geneBin_dict(args.annot_file,args.resolution)
	[highW0_bins, highW1_bins, highW2_bins] = get_topK_bins(Wval_list, K=args.topK_basis)

	highW_genes = [get_genelist(highW0_bins,gene_binDict), get_genelist(highW1_bins,gene_binDict), get_genelist(highW2_bins,gene_binDict)]
	genecountlist = [len(i) for i in highW_genes]	
	S_basis = np.argmax(genecountlist)
	print("[5/{}] Gene density measured".format(total_step))
	print("Gene density of top {} basis: {}".format(args.topK_basis, genecountlist))
	print("S-basis: {}".format(S_basis))

	outfilename = '{}/{}.{}.{}.W_basis.txt'.format(args.output_dir, args.header, args.resolution, S_basis)
	write_basis_file(outfilename, Wval_list, genecountlist)

	outfilename = '{}/{}.{}.Sbasis_genes.txt'.format(args.output_dir, args.header, args.resolution)
	write_gene_file(outfilename, highW_genes[S_basis])
	print("Job Done\n")

##
