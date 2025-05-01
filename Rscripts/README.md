# R scripts

`AlleleFrequency.R` processes data related to allele frequency. It utilises files named like `2612611001178601127QTL_V2_10k.csv` as input
It can plot number of QTLs, fixed QTLs, all allele frequencies, and changes in allele frequency per generation

`He+PhenotypePlotter.R` processes data related to the heritability and mean phenotype value. It utilises files like `2612611001178601127QTL_V2_He+Pheno_10k.csv` as input. 
It plots how heritability and allele frequency change through time. 

`LinkageDisequillibrium.R` calculates Linkage disequilibrium D, and D' (D prime). It utilised files like `QTL_V2_HaplotypeCount3528387796.csv`. One such file can be found in `genArchSelect/QTL_prototype/output_files/QTL_V2`. It contains information about which individuals contain which mutations. Go to `genArchSelect/QTL_prototype/scripts/QTL_V2.slim` to see how such a file can be output from the simulation.

`EffectSizevsDeltaP.R` plots the change in allele frequency per generation against the effect size. This is to test if higher effect sizes, change in frequency more.

The rest are files which were used to tweak and eventually arrive at these final files. `HeritabilityPlotter.R` and `LD_calc.R` contain most of this rough work. It maybe of importance to note that there is a function in `HeritabilityPlotter.R` called `SurcPlotter`, which does all the plotting for three graphs together, collates them and outputs it. 

`LD_calc_QTLs.R` is a script that calculates the LD with QTLs from the same and different chromosomes. Needs to be used with files from `~/QTLFixed`