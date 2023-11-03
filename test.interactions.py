import argparse
import sys
import os
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import zscore
from sklearn.cluster import KMeans
import seaborn as sns
import pickle

def read_spacemake(n):
	f = open(n)
	m = {}
	barcodes = set([])
	genes = set([])
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		t_key = (ll[0], ll[1])
		m[t_key] = int(ll[2])
		genes.add(ll[0])
		barcodes.add(ll[1])
	f.close()
	return m, genes, barcodes
def read_coord(n):
	f = open(n)
	Xcen = []
	cells = []
	for l in f:
		l = l.rstrip("\n")
		ll = l.split(",")
		cells.append(ll[0])
		Xcen.append((float(ll[-2]), float(ll[-1])))
	f.close()
	Xcen2 = np.empty((len(Xcen), 2), dtype="float32")
	for ind,(i,j) in enumerate(Xcen):
		Xcen2[ind, :] = [-1.0*i,j]
	Xcen = Xcen2
	return Xcen, cells

def read_expression(n, do_zscore=False):
	f = open(n)
	h = f.readline().rstrip("\n").split()
	h = [xh.replace(".", "-") for xh in h]
	num_gene = 0
	for l in f:
		l = l.rstrip("\n")
		num_gene+=1
	f.close()
	mat = np.empty((num_gene, len(h)), dtype="float32")
	f = open(n)
	f.readline()
	ig = 0
	genes = []
	for l in f:
		l = l.rstrip("\n")
		ll = l.split()
		gene = ll[0]
		values = [float(v) for v in ll[1:]]
		mat[ig,:] = values
		genes.append(gene)
		ig+=1
	f.close()
	if do_zscore:
		mat = zscore(mat, axis=1) #per cell
		#mat = zscore(mat, axis=0) #per gene
		
	return mat, h, genes

def read_genes(n):
	m = []
	f = open(n)
	for l in f:
		l = l.rstrip("\n")
		m.append(l)
	f.close()
	return m

def read_interactions(n):
	f = open(n)
	m = []
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		g1 = ll[2]
		g2 = ll[3]
		m.append((g1,g2))
	f.close()
	return m

def read_orthology(n):
	f = open(n)
	mouse_to_human, human_to_mouse = {}, {}
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		t_mouse = ll[0]
		t_human = ll[1]
		mouse_to_human.setdefault(t_mouse, [])
		mouse_to_human[t_mouse].append(t_human)
		human_to_mouse.setdefault(t_human, [])
		human_to_mouse[t_human].append(t_mouse)
	f.close()
	return mouse_to_human, human_to_mouse


if __name__=="__main__":
	parser = argparse.ArgumentParser(description="Prepare Gene Expression", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--mouse", dest="mouse", type=str, required=True, help="Xenomake mouse gene expression")
	parser.add_argument("--human", dest="human", type=str, required=True, help="Xenomake human gene expression")

	parser.add_argument("--barcode", dest="barcode", type=str, required=True, help="filtered barcodes and positions")
	parser.add_argument("--outdir", dest="outdir", type=str, required=True, help="output directory")
	args = parser.parse_args()


	m, genes, barcodes = read_spacemake(args.mouse) #m.keys(): gene,barcode "xenomake.filtered.mouse.txt"
	Xcen, Xcells = read_coord(args.barcode) #spatial_barcodes_visium_filtered.cvs

	map_cell = {}
	for ind,v in enumerate(Xcells):
		map_cell[v] = ind
	map_gene = {}
	for ind,v in enumerate(genes):
		map_gene[v] = ind

	#mouse matrix ==============================================
	mat = np.zeros((len(genes), len(Xcells)), dtype="float32")
	for t_gene,t_barcode in m:
		if t_barcode in map_cell and t_gene in map_gene:
			mat[map_gene[t_gene], map_cell[t_barcode]] = m[(t_gene, t_barcode)]

	m2, genes2, barcodes2 = read_spacemake(args.human) #m.keys(): gene,barcode #xenomake.filtered.human.txt
	map_gene2 = {}
	for ind,v in enumerate(genes2):
		map_gene2[v] = ind

	#human matrix =================================================
	mat2 = np.zeros((len(genes2), len(Xcells)), dtype="float32")
	for t_gene,t_barcode in m2:
		if t_barcode in map_cell and t_gene in map_gene2:
			mat2[map_gene2[t_gene], map_cell[t_barcode]] = m2[(t_gene, t_barcode)]

	if not os.path.exists(args.outdir):
		print("Creating output directory...")
		os.makedirs(args.outdir)
	fw = open("%s/map_genes.pkl" % args.outdir, "wb")
	pickle.dump([map_cell, map_gene, map_gene2, mat, mat2], fw)
	fw.close()
