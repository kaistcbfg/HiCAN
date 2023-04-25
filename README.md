# HiCAN
## Introduction
To systematically identify genomic regions establishing inter-chromosomal interactions from Hi-C contact maps, we devised a new computationalmethod named Hi-C interchromosomal contact map analysis with NMF (HiCAN).
The core design principle of HiCAN is mainly composed of three steps: construction of an intrachromosomal interaction-filled inter-chromosomal Hi-C contact map, projection of the contact map into three low-dimensional spaces with NMF, and annotation of S (speckle-associated)-, N (nucleolus-associated)-, and U (undefined)-basis based on gene density.
The numeric value of each entry (or, equivalently, each genomic locus) in an NMF basis indicates the degree to which the genomic locus belongs to the corresponding basis.
When we decomposed Hi-C contact maps with factorization rank 3, S and N bases values showed highly skewed distributions to specific genomic regions exclusive to each other,representing distinct modes of inter-chromosomal organizations. 
In contrast, U basis values followed a normal distribution across the entire genome, indicating genome-wide background signals. 
Importantly, the non-orthogonality and non-negativity of NMF bases enable effective deconvolution of the complex intermingled inter-chromosomal interactions into biological interpretable low-dimensional structures.
Thus, using HiCAN, we anticipated identifying three different modes of inter-chromosomal organizations.

## Installation
HiCAN requires following packages:
+ Numpy
+ Pandas
+ Matplotlib (optional, for visualization)
+ scikit-learn (optional, for Python NMF use)
+ rpy2 (optional, for R NMF use)
+ iced (optional, to run ICE normalization in script, https://github.com/hiclib/iced)

To use rpy2, 'NMF' package installed R environment settings are required.

## Usage
```
HiCAN: inter-chromosomal hub interaction caller

optional arguments:
  -h, --help            show this help message and exit
  --header HEADER       save file header
  --input-file INPUT_FILE
                        input *.gz file (required)
  --input-format INPUT_FORMAT
                        default pickle, format: text or pickle (numpy array pickle)
  --gzip-flag GZIP_FLAG
                        default True, *.gz input
  --resolution RESOLUTION
                        bin resolution (default 500kb
  --nmf-solver NMF_SOLVER
                        default R, format: R or Python
  --randstate-py RANDSTATE_PY
                        random state value, default None
  --output-dir OUTPUT_DIR
                        default ./output
  --genome-version GENOME_VERSION
                        default hg19, [hg19, hg38, mm10] available
  --genome-bindir GENOME_BINDIR
                        default ./genome_Bin. Look [genome].[chrname].[resolution].bin file in dir
  --chrname-list CHRNAME_LIST
                        default ./CHRLIST.txt List of chromosomes
  --coverage-cutoff COVERAGE_CUTOFF
                        default 1/3 of mean coverage (-1), should be > 0, row/col > thresh will remain
  --blacklist-file BLACKLIST_FILE
                        default False, remove corresponding rows/cols from input matrix
  --ICE-flag ICE_FLAG   default False, apply default ICE norm to matrix
  --annot-file ANNOT_FILE
                        gene annotation file
  --topK-basis TOPK_BASIS
                        default 500, get top K n-th basis to calc gene density
  --visualize-flag VISUALIZE_FLAG
                        Save contact maps to PDFs
  --visualize-vmax VISUALIZE_VMAX
                        default 50, visualize vmax value
  --visualize-output VISUALIZE_OUTPUT
                        PDF file save path
```

## Sample Data and Input Format

### Text foramt
Tab-separated Hi-C interaction table where 1st row (index 0) and column contains Hi-C bin coordinates.
Official example is in this format. Use following command to download.
```sh
wget https://dl.dropbox.com/s/ic6mmpyi5ozp3ns/IMR90_HiCAN_example.txt.gz
```
the example file is IMR90 cell line Hi-C data (hg38, chr1-22, 500kb resolution). It was processed by the blacklist removal, coverage filter, and ice normalization. Intra-chromosomal interactions were also removed (NaNs in table).
![input](https://dl.dropbox.com/s/u7hya3hyu881nwh/input_data_format.png)

Note that When using text input, the index information of row/col in text file is used and external index information given by --genome-version, --genome-bindir, and --chrname-list options are ignored.

### Numpy array

Pickled numpy array can be used as an input format.
The pickle file must contain only numpy arrays without index information.
External index information given by --genome-version, --genome-bindir, and --chrname-list options are used.

## Running Examples
### Run

Example HiCAN run script:
```bash
python HiCAN.py --header IMR90 --input-file IMR90_HiCAN_example.txt.gz  --input-format text --gzip-flag True --resolution 500000  --nmf-solver Python --visualize-vmax 20 --coverage-cutoff 0 --annot-file ./annot/gencode.v27.annotation_gene_protein_coding.gtf --topK-basis 500 --blacklist-file ./annot/hg38.blacklist.bed
```

Running following script will generate this output:
```
chrlist: ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']
Total 4818 bins loaded | resolution: 500000
row/cols from input file imported. --genome-bindir, --genome-version, and --chrname-list options ignored
IMR90 processing started....
[1/5] Data imported
Total 575 blacklist bins in resolution of 500000
ICE normalization OFF
[2/5] Input preprocessed. 4818 rows/cols remain after filtered
[3/5] Python NMF decomposition applied
[4/5] [1011, 1424, 2383] bins in each NMF basis
[5/5] Gene density measured
Gene density of top 500 basis: [1708, 5098, 1697]
S-basis: 1
Job Done
```

+ Note that Python NMF is faster, but results in our publication was originally generated with R NMFs.
+ Check options carefully to avoid misuse due to default settings
+ R NMF may generate extra logs 

### Output

Running the example command will generate two txt files and three PDFs in ./output and ./PDFs (can be changed by options).
+ ./output/IMR90.500000.1.W_basis.txt
+ ./output/IMR90.500000.Sbasis_genes.txt

W_basis.txt file contains three basis values of all valid bins.
1st row shows the gene density (genes are listed in Sbasis_genes.txt file) of each basis.
Index of S-basis (highest gene density) is marked in file name (1).

+ ./PDFs/IMR90_rawMat.pdf
+ ./PDFs/IMR90_filteredMat.pdf
+ ./PDFs/IMR90_reorderedMat.pdf

Reordered map shows re-ordered contact map according to the bins highest basis value.
Three well distinguished chunks are obsrevable.
![reorder](https://dl.dropbox.com/s/752qpejdqu6x6v1/IMR90_example_reordered.png)

As example file is already filtered, the rawMat and filteredMat is identical. 
Following figure shows the effect of preprocessing with another sample.
![filter](https://dl.dropbox.com/s/bggclyq9yzeryn5/filtering.png)


## Publication and Citation
> Jaegeon Joo, Sunghyun Cho, Sukbum Hong, Sunwoo Min, Kyukwang Kim, Rajeev Kumar, Jeong-Mo Choi, Yongdae Shin, and Inkyung Jung
> Probabilistic establishment of speckle-associated inter-chromosomal interactions
> *Nucleic Acids Research*, 2023;, gkad211
> doi: https://doi.org/10.1093/nar/gkad211

## Legacy
Early version of HiCAN source is available at:
http://junglab.kaist.ac.kr/Dataset/HiCAN_source.tar.gz

## License
For commercial use of the software, please contact the authors.