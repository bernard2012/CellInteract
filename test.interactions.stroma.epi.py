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
import argparse
import warnings

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
	parser = argparse.ArgumentParser(description="Binning Gene Expression (for Cross-Compartment interaction)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--inputdir", dest="inputdir", type=str, required=True, help="input directory which is outdir of previous step")
	parser.add_argument("--orthology", dest="orthology", type=str, required=True, help="human:mouse orthology file (get from github)")
	parser.add_argument("--interactlist", dest="interlist", type=str, required=True, help="ligand-receptor interaction list (get from github)")

	parser.add_argument("--nstart", dest="nstart", type=int, required=False, default=10, help="Kmeans nstart, recommend 1000 for accuracy or 10 for speed") #recommend nstart=1000 for most accuracy, nstart=10 for speed

	args = parser.parse_args()
	
	f = open("%s/map_genes.pkl" % args.inputdir, "rb")
	[map_cell, map_gene, map_gene2, mat, mat2] = pickle.load(f)
	f.close()

	#fw = open("map_genes.pkl", "wb")
	#pickle.dump([map_cell, map_gene, map_gene2, mat, mat2], fw)
	#fw.close()

	m2h, h2m = read_orthology(args.orthology) #human.mouse.orthologs.valid.txt 
	db = read_interactions(args.interlist) #cell.cell.interactions.good

	to_do_mouse = []
	to_do_human = []
	for g1,g2 in db:
		g1s = g1.split(";") #human gene
		g2s = g2.split(";") #human gene
		
		m1s, m2s = [], []
		for x in g1s:
			if x in h2m:
				for tt in h2m[x]:
					m1s.append(tt)
		for x in g2s:
			if x in h2m:
				for tt in h2m[x]:
					m2s.append(tt)

		m1s = sorted(list(set(m1s)))
		m2s = sorted(list(set(m2s)))

		if len(g1s)>0 and len(m2s)>0:
			for n_g1 in g1s:
				if n_g1 not in map_gene2:
					continue
				to_do_human.append(n_g1)
				for n_g2 in m2s:
					if n_g2 not in map_gene:
						continue
					to_do_mouse.append(n_g2)

		if len(m1s)>0 and len(g2s)>0:
			for n_g1 in m1s:
				if n_g1 not in map_gene:
					continue
				to_do_mouse.append(n_g1)
				for n_g2 in g2s:
					if n_g2 not in map_gene2:
						continue
					to_do_human.append(n_g2)
		

	to_do_mouse = sorted(list(set(to_do_mouse)))
	to_do_human = sorted(list(set(to_do_human)))
	print("mouse", len(to_do_mouse))
	print("human", len(to_do_human))
	member_mouse = np.empty((len(to_do_mouse), mat.shape[1]), dtype="int32")
	member_human = np.empty((len(to_do_human), mat2.shape[1]), dtype="int32")

	for ig,gg in enumerate(to_do_mouse):
		a_g1 = mat[map_gene[gg],:]
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			kmeans = KMeans(n_clusters=5, n_init=args.nstart).fit(np.transpose([a_g1]))
		member_mouse[ig,:] = kmeans.labels_
		print(ig, gg, kmeans.cluster_centers_)

	for ig,gg in enumerate(to_do_human):
		a_g1 = mat2[map_gene2[gg],:]
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			kmeans = KMeans(n_clusters=5, n_init=args.nstart).fit(np.transpose([a_g1]))
		member_human[ig,:] = kmeans.labels_
		print(ig, gg, kmeans.cluster_centers_)

	fw = open("%s/mouse_genes_stroma_epi_interactions_clusters.pkl" % args.inputdir, "wb")
	pickle.dump([to_do_mouse, member_mouse], fw)
	fw.close()

	fw = open("%s/human_genes_stroma_epi_interactions_clusters.pkl" % args.inputdir, "wb")
	pickle.dump([to_do_human, member_human], fw)
	fw.close()

