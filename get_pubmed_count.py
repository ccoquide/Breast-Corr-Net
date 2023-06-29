import os
from Bio import Entrez

# query function from the tutorial, worked fine on the first trial
def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='20',
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results

if __name__ == '__main__':

    # Get gene list into a list
    # Downloaded from :
    # https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
    # The protein_coding.txt file is generated with :
    # grep "^#\|protein-coding" Homo_sapiens.gene_info > protein_coding.txt
    
    all_genes = list()

    with open('C:\\Users\\alexo\\Documents\\JP_JL_CC\\protein_coding.txt') as filein:
        next(filein)
        for gene in filein:
            all_genes.append(gene.split('\t')[2])

    # get previously fetched pubmed count
    gene_already_done = list()
    tmp = open("C:\\Users\\alexo\\Documents\\JP_JL_CC\\gene_pubmed_count.txt", "r")
    for line in tmp:
        gene_already_done.append(line.split('\t')[0])
    tmp.close()

    # fetch and append into a file                           
    outputFile = open("C:\\Users\\alexo\\Documents\\JP_JL\\gene_pubmed_count.txt", "a")        
    index = 0
    for i in all_genes:
        if i not in gene_already_done:
            outputFile.write(i + '\t' + search(i)['Count'] + '\n')

        index += 1
        if index % 100 == 0 :
            print(round((index/len(all_genes))*100,2)+'%', end='\r')
    
    outputFile.close()
