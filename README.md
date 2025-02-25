# FuNTB

*`FuN-TB`* is a free, open-source Python tool designed for analyzing output from MTBSeq v.1.0. This standalone command-line tool enables the comparison of Single Nucleotide Polymorphisms (SNPs) across phenotypically diverse sets of Mycobacterium tuberculosis samples. By structuring the data into phenotype-centered networks, FUN-TB allows users to visualize which genes are altered and how these alterations correlate with specific phenotypes. The tool provides insights into gene alterations by representing them as nodes, with node size and edge width reflecting the strength of their relationship to particular phenotypes. 

*`FuN-TB`* consists of three core scripts:

1. *Variation Dictionary Creation*: Generates a dictionary of variations found in the `MTBSeq` output.
2. *Phenotype-Based Sample List Generation*: Generates a dictionary of variations found in the `MTBSeq` output.
3. *Phenotype-Centric and Gene-surrounded Network Structuring*: Constructs networks based on the relationship between phenotypic traits and gene alterations, making it easy to identify key genetic markers associated with specific phenotypes.
