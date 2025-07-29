# Bioinformatics - DADA2 🐠

### DADA2 pipeline

DADA2 is an R-based pipeline designed for high-resolution processing of amplicon sequencing data using **Amplicon Sequence Variants (ASVs)** insted of clustering into OTUs. 

## Workflow
1. **Filter and Trim** - Corrects sequencing errors to recover true biological sequences.
2. **Learn Error Rates** - Models sequencing errors from the data.
3. **Dereplication** - Combines identical reads for efficiency.
4. **Denoising** - Corrects errors and infers true sequences.
5. **Chimera Removal** - Filters out artificial recombinants.
6. **Taxonomy Assignment** - Classifies ASVs (e.g., with SILVA, Greengenes).

![image](https://github.com/user-attachments/assets/fcaa891c-dc5e-4977-b82c-25c49b19cdf4)

[DADA2 Documentation](https://benjjneb.github.io/dada2/tutorial.html)

