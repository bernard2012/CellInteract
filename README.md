## Running Cell-Cell Interaction Analysis

### Pre-requisites:

#### Input files:
- `human.mouse.orthologs.valid.txt` - human:mouse ortholog list (get from this repo)
- `cell.cell.interactions.good` - a list of ligand-receptor pairs (get from this repo)
- `xenomake.filtered.mouse.txt` - Xenomake mouse gene expression. [Example file](https://github.com/bernard2012/CellInteract.data/raw/main/xenomake.filtered.mouse.txt)
  - Gene, Barcode, UMI (positive) *tab-delimited*:
```
2900026A02Rik   AGGCCTGAGAATCTCG    2
2900026A02Rik   GTTATAATACGGTGAA    1
2900026A02Rik   ATCTCGTGAGCGAAAC    1
2900026A02Rik   CGCAATTACTTTCGGT    1
2900026A02Rik   ATTACCACACTGCCTG    2
2900026A02Rik   GACGCTTGCTTCTAAA    1
2900026A02Rik   ACGCTAGTGATACACT    2
2900026A02Rik   CTAGGGATAGGGACAA    1
2900026A02Rik   ACCATCCGCCAACTAG    1
2900026A02Rik   CGAAGTTGCTCTGTGT    1
2900026A02Rik   TTCCCGGCGCCAATAG    1
...
```
- `xenomake.filtered.human.txt` - Xenomake human gene expression [Example file](https://github.com/bernard2012/CellInteract.data/raw/main/xenomake.filtered.human.txt)
  - Gene, Barcode, UMI (positive) *tab-delimited*
```
MKLN1   AACTTGCCCGTATGCA    1
MKLN1   GTTGGTCATGCTATCC    1
MKLN1   ATTTGCCTAGTTACGA    1
MKLN1   CACTGACGATTGTGGA    1
MKLN1   GATCGACACTATCTGA    1
MKLN1   TTGACATGAACGTGGA    2
MKLN1   CTAGTAGAAAGGGATT    1
MKLN1   GGGTCAGGAGCTAGAT    1
MKLN1   AAAGTGTGATTTATCT    1
MKLN1   TAACGCTTTGAGAGCG    1
MKLN1   TTGCCATAGCCCGCTC    2
MKLN1   GCAGCACACAGCCCAG    1
MKLN1   GTATAGGACTCAGTAG    1
MKLN1   GATCCCTTTATACTGC    2
...
```
- `spatial_barcodes_visium_filtered.csv` - a list of Visium filtered barcodes (i.e. in-tissue barcodes) and their positions. [Example file](https://github.com/bernard2012/CellInteract.data/raw/main/spatial_barcodes_visium_filtered.csv)
  - barcode,x,y
```
AAGAGCTCTTTATCGG,50,68.92307692307692
AACGGCCATCTCCGGT,88,49.230769230769226
ATACTACCCGTACCAC,22,29.538461538461533
GTGGCGGTCCCAGCGT,107,60.71794871794872
TTACAGACCTAAATGA,31,83.6923076923077
GCGCCGTTCCACGATA,62,59.07692307692308
CAAGATATTATAACGT,38,52.51282051282051
GAGGCCCGACTCCGCA,104,36.1025641025641
ATGTGAAAGCCTAATG,105,50.87179487179488
GACCGACTGAAGCGTC,105,70.56410256410257
```

#### Step 1:
Prepare gene expression matrix. 
Running `test.interactions.py`:
```
python3 test.interactions.py --mouse xenomake.filtered.mouse.txt --human xenomake.filtered.human.txt --barcode spatial_barcodes_visium_filtered.csv --outdir out
```

#### Step 2:

Discretize gene expression of each gene into 5 bins by K-means: 0 (lowest expression), 1, 2, 3, 4 (highest expression)

- Running `test.interactions.step2.py` on mouse
  - Output file:
`mouse_genes_cell_interactions_clusters.pkl` in `--inputdir out`
```
python3.6 test.interactions.step2.py --inputdir out --orthology human.mouse.orthologs.valid.txt --interactlist cell.cell.interactions.good --nstart 10 --do-organism mouse
```
- Running `test.interactions.step2.py` on human
  - Output file:
`human_genes_cell_interactions_clusters.pkl` in `--inputdir out`
```
python3.6 test.interactions.step2.py --inputdir out --orthology human.mouse.orthologs.valid.txt --interactlist cell.cell.interactions.good --nstart 10 --do-organism human
```

#### Step 3:

Run cell-cell interactions using Giotto's algorithm. Returns significant ligand-receptor pairs.

- Run `calc.distance.giotto.py` to get **epithelial-epithelial** interactions (setting `--do-organism human`)
```
python3 calc.distance.giotto.py --inputdir out --barcode spatial_barcodes_visium_filtered.csv --orthology human.mouse.orthologs.valid.txt --interactlist cell.cell.interactions.good --do-organism human
```
- Run `calc.distance.giotto.py` to get **stromal-stromal** interactions (setting `--do-organism mouse`)
```
python3 calc.distance.giotto.py --inputdir out --barcode spatial_barcodes_visium_filtered.csv --orthology human.mouse.orthologs.valid.txt --interactlist cell.cell.interactions.good --do-organism mouse
```
- Output files are located in the same directory as `--inputdir out`.
- Description of output file `mouse_stroma_stroma_interactions.txt`:
```
Tnc Itgb1   33.802794
Fn1 Itgb1   32.621910
Fn1 Itga5   27.240370
Col6a1  Itgb1   27.189511
Spp1    Itgb1   26.994819
Col1a2  Itgb1   26.579064
Fbn1    Itgb1   26.182035
Col3a1  Itgb1   25.751852
Col16a1 Itgb1   25.467074
Col8a1  Itgb1   25.029367
Apoe    Tyrobp  24.560072
Apoe    Trem2   23.619534
Col5a1  Itgb1   23.100303
Col5a2  Itgb1   23.098580
Thbs1   Itgb1   22.394043
Col6a3  Itgb1   21.719845
Tnc Itgav   21.398422
...
```
- **First column**: ligand; 
- **second column**: receptor; 
- **third column**: interaction z-score.


#### Step 4:

**Cross-compartment** (i.e. stromal-epithelial and epithelial-stromal) interactions.

Binning gene expression for cross-compartment interactions
```
python3.6 test.interactions.stroma.epi.py --inputdir out --orthology human.mouse.orthologs.valid.txt --interactlist cell.cell.interactions.good --nstart 10
```

Test ligand receptor interactions:
```
python3.6 calc.distance.stroma.epi.giotto.py --inputdir out --barcode ../spatial_barcodes_visium_filtered.csv --orthology human.mouse.orthologs.valid.txt --interactlist cell.cell.interactions.good
```

- Output files are located in directory specified in `--inputdir out`. `epithelium_stroma_interactions.txt` and `stroma_epithelium_interactions.txt`.
- Description of output (`epithelium_stroma_interactions.txt`):
```
COL1A1  Itgb1   13.333940
COL1A1  Itga11  8.151005
COL1A1  Itga1   6.770442
PTGES3  Ptger3  5.002067
COL8A1  Itga11  3.195053
WNT2B   Sfrp4   2.929544
TNC Itgb3   2.905398
ARTN    Ret 2.858487
FN1 Itgb3   2.787136
COL4A1  Itga2   2.720209
CDH1    Itgb7   2.710135
CADM1   Cadm1   2.689628
SLC6A12 Gabbr1  2.642121
RSPO1   Kremen1 2.591782
PTGES2  Ptger3  2.571076
TNC Itga9   2.526410
JAG2    Notch4  2.470840
BMP2    Acvr2a  2.464271
WNT3    Fzd5    2.319953
TF  Tfrc    2.239750
HSD17B12    Ar  2.159742
RARRES2 Ccrl2   2.151600
DHCR24  Rora    2.014379
COL8A1  Itgb1   1.999701
WNT4    Fzd4    1.997669
ICAM1   Itgal   1.927682
TSLP    Crlf2   1.925742
TSLP    Il7r    1.830379
...
```
- **First column**: ligand (from human)
- **second column**: receptor (from mouse); 
- **third column**: interaction z-score.

#### Output file list
- `epithelium_stroma_interactions.txt` [Example output](https://media.githubusercontent.com/media/bernard2012/CellInteract.data/main/example.out/epithelium_stroma_interactions.txt)
- `stroma_epithelium_interactions.txt` [Example output](https://media.githubusercontent.com/media/bernard2012/CellInteract.data/main/example.out/stroma_epithelium_interactions.txt)
- `mouse_stroma_stroma_interactions.txt` [Example output](https://media.githubusercontent.com/media/bernard2012/CellInteract.data/main/example.out/mouse_stroma_stroma_interactions.txt)
- `human_epithelium_epithelium_interactions.txt` [Example output](https://media.githubusercontent.com/media/bernard2012/CellInteract.data/main/example.out/human_epithelium_epithelium_interactions.txt)
