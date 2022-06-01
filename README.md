# Breast-Corr-Net

# Networks files, nodes' information and preliminary results associated to Methylated-DNA X RNA pearson correlations in the context of about 800 tissues from breast cancer.

## Dependencies and Setup

1. Cytoscape version >= 3.8.2 (https://cytoscape.org/)

- Files `CC_col/RNA.cys` and `CC_row/MethDNA.cys` has to be opened with Cytoscape and contains network style presets for RNA-RNA and DNA-DNA correlations sub-networks visualisation respectively.

## Filenames information

### Files required for network's visualization using Cytoscape

1. `CorMethRna-0.6-trans-proof-proj-[col/row]-[type]-new-CC-[n_cc].net` : List of links related to the `[n_cc]`-th cluster extracted from the projection over the col (RNA) or row (Methylated DNA) space containing either only positive (`[type]` = positive) or negative (`[type]` = negative) complementary correlations or both of them (`[type]` = mixed).
2. `CorMethRna-0.6-trans-proof-proj-[col/row]-[type]-new-CC-[n_cc]-notinproj.net` : List of links between nodes present in `[n_cc]`-th cluster of type `[type]` associated to the projection over col/row space  that are only present in the observed network (the one constructed from RNA-RNA or DNA-DNA observed correlations).
3. `CorMethRna-0.6-trans-proof-proj-[col/row]-[type]-new-CC-[n_cc]-notObserved.net`  : List of links not present in the observed network

These files are formated as following (with first line containing Column names)
	
        C1 C2 C3 C4
        S1 T1 W1 IP1
        S2 T2 W2 IP2
        S3 T3 W3 IP3
        ...
	
Where `Si` and `Ti` are the nodes inedex (integer) of the i-th link's source and target nodes, `Wi` its weight (integer for almost all `[*].net` files, exepted for `[*]-notinproj.net` where it is a float) and `IPi` a boolean (1 if the link is present in the observed network, 0 else).

4. `CorMethRna-0.6-trans-proof-proj-[col/row]-[type]-new-CC-[n_cc].nodes` : List of network's nodes information such as key ID (integerer), node name (ex: `ENSG00000054598`). This nodes file can be used with the `[*].net`, `[*]-notinproj.net` and `[*]-notObserved.net` network files.

The nodes files are formated as follows
	
       C1 C2 C3 C4 C5 ...
       ID1 NAME1 IP1 EXIF11 EXIF21 ...
       ID2 NAME2 IP2 EXIF12 EXIF22 ...
       ID3 NAME3 IP3 EXIF13 EXIF23 ...
       ...
	
Where `IDi` is the node index (integer) used in the associated `[*].net` file, `NAMEi` is the corresponding node's name, `IPi` is a boolean (1 if the node is not present in the observed network, 0 else), and finally `EXIF1i` is the first extra option such as the chromosome name, or other information related to either RNA or Methylated-DNA.


**Note:** An empty `[*]-notinproj.net` or `[*]-notObserved.net` file means that such correlations don't exist.
### Other files

1. `CorMethRna-0.6-trans-proof-proj-[col/row]-[type]-new-CC-[n_cc].info` : Information of the network of the same `*.net` name. Lists the percentage of nodes being present in the observed network (the one constructed from RNA-RNA or DNA-DNA observed correlations) + the percentage of links being present in the observed network.
2. `[col/row]-[type]-netdensity-sorted.info` : Information of the different clusters related to the projection over col/row space containing correlations of type `[type]` (ex: Linkage Density, Number of nodes ...)
3. `Row-MethDNA.txt`: List of all nodes related to the row space (Methylated DNA) and related information.
4. `Col-RNA.txt`: List  of all nodes related to the col space (RNA) and related information.
5. `report.pdf` PDF file containing preliminary results and a draft of Material and Methods.
5. `sparse-[*]-net-example` and `small-[*]-net-example` are text files with information of the cluster of interests that are presented in the `report.pdf` file

### Directories

1. `CC_col`:  Contains all clusters related to the `col` space projection.
2. `CC_row`: Contains all clusters related to the `row` space projection.
3. `POSITIVE`: Sub-directory from `CC_col` or `CC_row`containing all clusters related to `[type]` = `positive` correlations.
4. `NEGATIVE`: Sub-directory from `CC_col` or `CC_row`containing all clusters related to `[type]` = `negative` correlations.
5. `MIX`: Sub-directory from `CC_col` or `CC_row`containing all clusters related to `[type]` = `mixed` correlations.
6. `REPORT`: Contains latex file and PDF of the report + vector images used as figures.
