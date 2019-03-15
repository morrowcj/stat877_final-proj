# data file descriptions

All files with the `gwa-files/wisasp_gwa-data` prefix are the data files. Below
are descriptions of each file according to the file extention

* `.bed`: binary file containing genotype information for each genet*SNP

* `.fam`: metadata to be used with .bed files w/ six columns:tags: NULL

    - Family ID: identical to the individual ID - there is no population structure, 
    individuals are thoroughly unrelated.
    
    - Individual ID: Genet identifier (identical to Family ID)
    
    - Paternal ID: NA
    
    - Maternal ID: NA
    
    - Sex: NA - sex provided in phenotype data
    
    - Phenotype: NA - provided in phenotype files
    
* `.map` and `.bim`: map file with SNP meta dat - contains 4, and 6 columns
respectively (2 extra for `.bim`).

    - chromosome: all = 1, chromosome position is unknown. scaffoldings were 
    used instead
    
    - SNP: the SNP identifier of the format <scaffoldID>:<bp position>
    
    - bp: base pair position on the 'chromosome'
    
    - allele1: name of the reference allele - `.bim` only
    
    - allele2: name of the non-reference allele - `.bim` only

* `.ped` (not contained in repo - too large): phenotype data for each genet at each SNP. The first 6 columns are the 
same as the `.fam` file and the remaining columns are the alleles for each snp
(i.e. snp1_allele1, snp1_allele2, snp2_allele1, snp2_allele2, ...)

* `.nosex`: file containing list of genets with no sex info (all of them)
