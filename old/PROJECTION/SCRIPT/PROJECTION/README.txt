The python script buildSignedProjNet.py takes in input two parameters, the first one is the name of the 
correlation matrix file, for instance ADN-ARN.dat, the second is an option.
There are three possible options:
	opt 1 : A float in the range [0,1], it will be considered as a threshold when the correlation matrix will be loaded.
	opt 2 : "brute-force", in this case when the projected ARN/ARN (ADN/ADN) correlations are calculated, only ADN-ARN correlations >= 1/sqrt(2) will be taken into account. By doing this we ensure those projected correlations to be realistic (see the ref [1] from the report)
	opt 2 : "trans-proof", in this case the projected correlation will be also realistic in the sens of correlation transitivity, but without setting any threshold such as it is the case with "brute-force" option. With "trans-proof" option make the transitivity test for each pair \rho_1 and \rho_2 of correlations.

The input file is formated as follows:
	[Number of column]	[Number of rows]
	[Number of non-zero entries]
	[Column index]	[Row index]	[Correlation]
	[Column index]	[Row index]	[Correlation]
	...
The two firs rows are related to the bipartite network information, number of column nodes, row nodes and links
All other rows describe each pair of interaction, first column index, then row index then the correlation coeff. Note that indices are in the range [1,N] where N is the number of row (col) nodes.
Finally, each element is separated by a tabulation

If you try with the file from data_test, ADN-ARN.dat, column nodes corresponds to ARN and row nodes to ADN.

To execute the python script use the command: python3.7  buildSignedProjNet.py ADN-ARN.dat [opt] (opt can be a positive float, or "brute-force", or "trans-proof"
Ex:  python3.7  buildSignedProjNet.py ADN-ARN.dat 0.4 

The output files consists in:
	For each projection we have :
		[filename].dat : Contains the projected network and its information
		[Total Number of ARN (ADN) in the dataset ]
		[Number of links in the projection]
		[Number of positive projected correlation]	[Number of negative projected correlation]
		[node index]	[node index]	[Positive projected correlation intensity]	[Negative projected correlation intensity]
		[node index]	[node index]	[Positive projected correlation intesity]	[Negative projected correlation intensity]
		...
The three first rows give information about the projected network, then all other rows consist in the pair of interacting nodes with the weight associated to each possible link sign (positive / negative). Note that node index is in the rangfe [1,N], where N is the total number of node in the global bipartite network (total number of ADN or ARN).

		[filename].info :
		[Total intensity of negative projected correlations]
		[Total intensity of positive projected correlations]
		[Total intensity of mixed sign projected correlations]
		[number of nodes interacting positively and negatively with other nodes]
		[number of links in the projection]
