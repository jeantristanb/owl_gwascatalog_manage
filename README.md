# extraction of phenotype from OWL and gwas catalog linkerd to cancer 
## short algorithm
* EFO extraction link to cancer : 
 * extraction from owl format of all phenotype descent of `http://purl.obolibrary.org/obo/MONDO_0045024` (categorie : cancer or benign tumor) from `datai/efo.owl.gz` using bin/clusteriseowl.py
   * file contains categorie : `data_change/MONDO_0045024.descent`
* GWAS catalog :
 * download of gwas catalog studies file and  ancestrality (`datai/gwas-catalog-v1.0.3-ancestries-r2022-07-09.tsv`, `gwas-catalog-v1.0.3-studies-r2022-07-09.tsv`):
 * used EFO categorie to select cancer link   using `bin/gwascat.r`
 * plot of distribution by ancestrality of individuals

