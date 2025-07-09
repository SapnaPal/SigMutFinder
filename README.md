# SigMutFinder
Significant Mutation Finder - A web server to map significant cancer mutations on kinases and PPIs 
Some of the code for this tool is hidden due to privacy issues. It will become public soon after publication. COMING SOON!!!☺️


# Implementation of the bioinformatics pipeline to provide mutation analysis for kinases and PPIs 
To systematically identify and characterize cancer-associated mutations within kinase domains and their protein-protein interaction (PPI) interfaces, we applied our newly developed in-house computational pipeline, SigMutFinder. The SigMutFinder maps mutations from the input MAF files onto relevant kinase genes and PPI interfaces. The application of SigMutFinder includes: 

## Mutation mapping on kinases and PPI interface

Initially, the pipeline maps somatic mutations from input MAF files onto relevant kinase genes and protein-protein interaction interfaces. Subsequently, it also integrates the pre- computed pan-cancer mutation data derived from the statistical analysis of our study into the input mutations. This process outputs frequently mutated kinase genes, high-frequency mutations in specific cancer types, enriched pathways among top kinase genes, and potential oncoPPIs with significantly mutated interface residues (A).



## De Novo Statistical Analysis of Input Cohorts:

Furthermore, when applied to our specific cohort of input MAF files (especially for a large number of samples), Kin_PPI_MutFinder can also perform a de novo statistical analysis. This enables the user to directly identify and validate the presence of frequently mutated kinase genes, highly frequent mutations in specific cancer types within our dataset, enriched pathways in top kinase genes, and novel oncoPPIs with significantly mutated interface residues among the analysed input MAF files (B). The typical examples of output files produced by Kin_PPI_MutFinder is depicted in Figure 12C.



