# Scripts

QTL_V2_2.slim is the code for simulating a nucleotide-based polygenic. 
This is the base code, which is adapted for different selection scenarios and to run burn ins. It outputs information into two files:
- He+Pheno: a csv file with information on heritability and mean phenotype through generations
- Allele: a csv file with information about what QTLs are present in which generation and in what frequency

More info on the output files can be found in `output_files`

QTL_V2_BurnIn.slim is the same as above but adapted with the following. There is no selection present. Once the population reached drift-mutation balance, which is when the `calcHeterozygosity();` of the population reaches theta (4NMu), the population state is saved. This is used to generate burnins.

The file `withSelection` contains variations of QTL_V2_2.slim, but loads in burnt-in population states and puts them under various selection treatments. Looking at the file will show which burnin was used. All burnins are present in `output_files` under the appropriate parameters.

