from requests import get
from re import findall
import gzip
from os import chdir, path, getcwd

chdir('C:\\Users\\alexo\\Documents\\JP_JL')
outf = 'gene_pubmed_count_2.txt'
# Create the list to fill
genes_names = []
genes_ID = []
genes_aliases = []
genes_counts = []

with gzip.open('Homo_sapiens.gene_info.gz', 'rb') as f:
    next(f)
    for line in f:
        line = line.decode().strip().split('\t')
        # remove untranscripted regions
        if line[9] in ['snoRNA','snRNA','ncRNA','protein-coding']:
            genes_ID.append(line[1])
            genes_names.append(line[2])
            genes_aliases.append(line[4])
            genes_counts.append('0')

# header
if not path.exists(outf):
    outputFile = open('gene_pubmed_count_2.txt', 'w')
    outputFile.write('ID\tgene name\taliases\tpubmed_count\n')
else:
    # Check the genes already counted
    # This has to be done since the process is >10h and unstable (web requests)
    tmp = open(outf, 'r')
    for line in tmp:
        line = line.strip().split('\t')
        if line[0] in genes_ID:
            genes_counts[genes_ID.index(line[0])] = line[3]
    tmp.close()

    outputFile = open('gene_pubmed_count_2.txt', 'a')
    
geneUrl = 'https://www.ncbi.nlm.nih.gov/gene/?term=' # require gene name
pubmedUrl = 'https://pubmed.ncbi.nlm.nih.gov/?show_snippets=off&size=10&linkname=gene_pubmed&from_uid=' # require gene ID



index = 0
for i in range(len(genes_ID)):
    if genes_counts[i] == '0':
        
        ## Get data from pubmed
        # res = get(pubmedUrl + genes_ID[i])
        # pubmedCount = findall('class=\"value\">[0-9]+,?[0-9]*', res.text)

        ## Get data from gene
        res = get(geneUrl + genes_names[i])
        pubmedCount = findall('pubmed-count-link\">\([0-9]+,?[0-9]*', res.text)
        
        # Check that the gene exist in pubmed
        if len(pubmedCount) > 0:
            pubmedCount = findall('[0-9]+,?[0-9]*', pubmedCount[0])[0]
            genes_counts[i] = pubmedCount.replace(',','')
            outputFile.write(genes_ID[i]+'\t'+genes_names[i]+'\t'+genes_aliases[i]+'\t'+genes_counts[i]+'\n')
        else:
            if len(findall('Found 1 result',res.text)) > 0:
                genes_counts[i] = '1'
                outputFile.write(genes_ID[i]+'\t'+genes_names[i]+'\t'+genes_aliases[i]+'\t'+genes_counts[i]+'\n')
            else:
                outputFile.write(genes_ID[i]+'\t'+genes_names[i]+'\t'+genes_aliases[i]+'\t'+genes_counts[i]+'\n')

        if index % 30 == 0 :
            print(genes_ID[i]+'\t'+genes_names[i]+'\t'+genes_aliases[i]+'\t'+genes_counts[i])
        
        if index % 100 == 0 :
            print(str(round(index*100/len(genes_ID),1))+'%')
    index +=1
outputFile.close()
##tmp = []'
##async def scrape(gene):
##    try:
##        res = requests.get("https://www.ncbi.nlm.nih.gov/gene/?term=" + gene)
##        pubmedCount = re.findall('pubmed-count.*\)', res.text)[0]
##        pubmedCount = pubmedCount.replace('pubmed-count-link\">(','').replace(')','').replace(',','')
##        tmp.append(gene + '\t' + pubmedCount)
##    except:
##        tmp.append(gene + '\t' + "pubmed not availble")
##
##async def main(): 
##	tasks = [ 
##		scrape(i) 
##		for i in all_genes
##	] 
##	await asyncio.gather(*tasks) 
##
##
##asyncio.run(main())
##
##outputFile.write('\n'.join(tmp))
##
##outputFile.close()
##
##
##        
##pat = list(range(0,len(all_genes),1000))
##pat.append(len(all_genes))
##
##bob = all_genes[1:10]
##
##for i in bob:
##    scrape(i)

