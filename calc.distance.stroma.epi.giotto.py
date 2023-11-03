import sys
import os
import re
from numpy import random
from scipy import stats
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
from scipy.stats import zscore
from scipy.spatial.distance import squareform, pdist
from scipy.stats import percentileofscore
from operator import itemgetter, attrgetter
from sklearn.cluster import KMeans
import seaborn as sns
import pickle
from scipy.spatial import Delaunay
import warnings
import argparse

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

def rank_transform_matrix(mat, rbp_p = 0.99, reverse=True):
	dim1 = mat.shape[0]
	dim2 = mat.shape[1]
	rank_forward = np.empty([dim1, dim2])

	print("Start ranking forward...")
	for c1 in range(dim1):
		rd = scipy.stats.rankdata(mat[c1,:])
		if reverse==True:
			rd = dim2 - rd + 1
		rank_forward[c1, :] = rd + 1
		if c1%1000==0:
			print("Done %d" % c1)
	print("Finished ranking forward...")
	rank_backward = np.empty([dim1, dim2])

	print("Start ranking backward...")
	for c1 in range(dim2):
		rd = scipy.stats.rankdata(mat[:,c1])
		if reverse==True:
			rd = dim1 - rd + 1
		rank_backward[:, c1] = rd + 1
		if c1%1000==0:
			print("Done %d" % c1)

	print("Finished ranking backward...")
	mutual_rank_rbp = np.empty([dim1, dim2])
	mutual_rank = np.empty([dim1, dim2])
	print("Calculate mutual rank...")
	ma = np.sqrt(np.multiply(rank_forward, rank_backward))
	print("Calculate exponential transform...")
	mutual_rank_rbp = np.multiply(1-rbp_p, np.power(rbp_p, np.subtract(ma, 1)))
	print("Finished exponential transform...")
	mutual_rank = ma

	dissimilarity = np.empty([dim1, dim2])
	print("Calculate dissimilarity...")
	dissimilarity = np.subtract(1, np.divide(mutual_rank_rbp, 1-rbp_p))
	print("Finished dissimilarity...")
	return dissimilarity


def analyze_member(member, expr):
	ranking = []
	max_per = 0.10
	min_spot = 10
	for i in range(5):
		m = np.where(member==i)[0]
		avg = np.mean(expr[m])
		if np.isnan(avg): continue
		if avg==0: continue
		ranking.append((i, avg, m.shape[0]))
	ranking = sorted(ranking, key=itemgetter(1), reverse=True)

	#use running_sum
	good_ranking = []
	tot = member.shape[0]
	running_sum = 0
	for i,j,k in ranking:
		tmp = running_sum + k
		if tmp/tot>max_per:
			break
		good_ranking.append((i, j, k, max_per - tmp/tot))
		running_sum += k

	#good_ranking = sorted(good_ranking, key=itemgetter(3))

	'''
	ind = good_ranking[0]
	tm = np.zeros(member.shape[0], dtype="int32")
	m = np.where(member==ind[0])[0]
	tm[m] = 1
	t_new = np.copy(tm)
	m2 = np.where(t_new==1)[0]
	if m2.shape[0]>=min_spot:
		good_ind.append(m2)
	'''
	
	good_ind = []
	tm = np.zeros(member.shape[0], dtype="int32")
	for i,j,k,l in good_ranking:
		m = np.where(member==i)[0]
		tm[m] = 1
		t_new = np.copy(tm)
		m2 = np.where(t_new==1)[0]
		if m2.shape[0]>=min_spot:
			good_ind.append((m2, l))
			#good_ind.append(m2)

	good_ind = sorted(good_ind, key=itemgetter(1))
	if len(good_ind)==0:
		return []

	#print(good_ind)
	#do not use running_sum
	'''
	good_ranking = []
	tot = member.shape[0]
	#running_sum = 0
	for i,j,k in ranking:
		#tmp = running_sum + k
		tmp = k
		if tmp/tot>max_per:
			break
		good_ranking.append((i,j,k))
		#running_sum += k
	good_ind = []
	tm = np.zeros(member.shape[0], dtype="int32")
	for i,j,k in good_ranking:
		m = np.where(member==i)[0]
		tm[m] = 1
		t_new = np.copy(tm)
		m2 = np.where(t_new==1)[0]
		if m2.shape[0]>=min_spot:
			good_ind.append(m2)
	'''
	#print(len(good_ind[0][0]))
	return [good_ind[x][0] for x in range(len(good_ind))]
	#return [good_ind[0][0]]
	#return good_ind

def count_edges(ind1, ind2, neighbors):
	edges = set([])
	for i1 in ind1:
		nn = set(ind2) & neighbors[i1]
		for x_nn in nn:
			edges.add(tuple(sorted([i1, x_nn])))
	#for i2 in ind2:
	#	nn = set(ind1) & neighbors[i2]
	#	for x_nn in nn:
	#		edges.add(tuple(sorted([i2, x_nn])))
	return len(edges)


if __name__=="__main__":
	parser = argparse.ArgumentParser(description="Test interactions", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--inputdir", dest="inputdir", type=str, required=True, help="input directory which is outdir of previous step")
	parser.add_argument("--barcode", dest="barcode", type=str, required=True, help="filtered barcodes and positions")
	parser.add_argument("--orthology", dest="orthology", type=str, required=True, help="human:mouse orthology file (get from github)")
	parser.add_argument("--interactlist", dest="interlist", type=str, required=True, help="ligand-receptor interaction list (get from github)")
	args = parser.parse_args()
	
	f = open("%s/map_genes.pkl" % args.inputdir, "rb")
	[map_cell, map_gene, map_gene2, mat, mat2] = pickle.load(f)
	f.close()
	Xcen, Xcells = read_coord(args.barcode) #spatial_barcodes_visium_filtered.csv

	tri = Delaunay(Xcen)
	print(tri.simplices.shape)

	neighbors = {}
	for xt in range(tri.simplices.shape[0]):
		points = tri.simplices[xt,:]
		neighbors.setdefault(points[0], set([]))
		neighbors[points[0]].add(points[1])
		neighbors[points[0]].add(points[2])
		neighbors.setdefault(points[1], set([]))
		neighbors[points[1]].add(points[0])
		neighbors[points[1]].add(points[2])
		neighbors.setdefault(points[2], set([]))
		neighbors[points[2]].add(points[0])
		neighbors[points[2]].add(points[1])

	adj = np.empty((Xcen.shape[0], Xcen.shape[0]), dtype="int32")
	for n in neighbors:
		for m in neighbors[n]:
			adj[n, m] = 1
			adj[m, n] = 1

	m2h, h2m = read_orthology(args.orthology) #human.mouse.orthologs.valid.txt
	db = read_interactions(args.interlist) #cell.cell.interactions.good

	f = open("%s/human_genes_stroma_epi_interactions_clusters.pkl" % args.inputdir, "rb")
	[to_do_human, member_human] = pickle.load(f)
	f.close()
	map_todo_human = {}
	for ig,g in enumerate(to_do_human):
		map_todo_human[g] = ig
	

	f = open("%s/mouse_genes_stroma_epi_interactions_clusters.pkl" % args.inputdir, "rb")
	[to_do_mouse, member_mouse] = pickle.load(f)
	f.close()
	map_todo_mouse = {}
	for ig,g in enumerate(to_do_mouse):
		map_todo_mouse[g] = ig


	tot_spot = member_human.shape[1]
	#to_do_mouse = []
	pairs_hm = []
	pairs_mh = []
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
				for n_g2 in m2s:
					if n_g2 not in map_gene:
						continue
					pairs_hm.append((n_g1, n_g2))

		if len(m1s)>0 and len(g2s)>0:
			for n_g1 in m1s:
				if n_g1 not in map_gene:
					continue
				for n_g2 in g2s:
					if n_g2 not in map_gene2:
						continue
					pairs_mh.append((n_g1, n_g2))

	pairs_hm = set(pairs_hm)
	pairs_mh = set(pairs_mh)

	outlines = []
	for n_g1,n_g2 in pairs_hm:
		m_g1 = member_human[map_todo_human[n_g1],:]
		m_g2 = member_mouse[map_todo_mouse[n_g2],:]
		a_g1 = mat2[map_gene2[n_g1],:]
		a_g2 = mat[map_gene[n_g2],:]
		indices1 = analyze_member(m_g1, a_g1)
		indices2 = analyze_member(m_g2, a_g2)
		if len(indices1)==0 or len(indices2)==0:
			continue
		vals = []
		for ind1 in indices1:
			for ind2 in indices2:
				#random
				rand_dist = []
				for xi in range(1000):
					rand1 = random.randint(tot_spot, size=(ind1.shape[0]))
					rand2 = random.randint(tot_spot, size=(ind2.shape[0]))
					rand_dist.append(np.sum(adj[np.ix_(rand1, rand2)]))
				across_dist = np.sum(adj[np.ix_(ind1, ind2)])
				t_mean = np.mean(rand_dist)
				t_std = np.std(rand_dist)
				vals.append((across_dist - t_mean)/t_std)
				#vals.append(stats.percentileofscore(rand_dist, across_dist))
		print(n_g1, n_g2, np.mean(vals))
		outlines.append((n_g1, n_g2, np.mean(vals)))

	outlines = sorted(outlines, key=itemgetter(2), reverse=True)
	fw =open("%s/epithelium_stroma_interactions.txt" % args.inputdir, "w")
	for i,j,k in outlines:
		fw.write("%s\t%s\t%f\n" % (i,j,k))
	fw.close()

	outlines = []
	for n_g1,n_g2 in pairs_mh:
		m_g1 = member_mouse[map_todo_mouse[n_g1],:]
		m_g2 = member_human[map_todo_human[n_g2],:]
		a_g1 = mat[map_gene[n_g1],:]
		a_g2 = mat2[map_gene2[n_g2],:]
		indices1 = analyze_member(m_g1, a_g1)
		indices2 = analyze_member(m_g2, a_g2)
		if len(indices1)==0 or len(indices2)==0:
			continue
		vals = []
		for ind1 in indices1:
			for ind2 in indices2:
				#random
				rand_dist = []
				for xi in range(1000):
					rand1 = random.randint(tot_spot, size=(ind1.shape[0]))
					rand2 = random.randint(tot_spot, size=(ind2.shape[0]))
					rand_dist.append(np.sum(adj[np.ix_(rand1, rand2)]))
				across_dist = np.sum(adj[np.ix_(ind1, ind2)])
				t_mean = np.mean(rand_dist)
				t_std = np.std(rand_dist)
				vals.append((across_dist - t_mean)/t_std)
				#vals.append(stats.percentileofscore(rand_dist, across_dist))
		#print(n_g1, n_g2, vals)
		print(n_g1, n_g2, np.mean(vals))
		outlines.append((n_g1, n_g2, np.mean(vals)))

	outlines = sorted(outlines, key=itemgetter(2), reverse=True)
	fw =open("%s/stroma_epithelium_interactions.txt" % args.inputdir, "w")
	for i,j,k in outlines:
		fw.write("%s\t%s\t%f\n" % (i,j,k))
	fw.close()



