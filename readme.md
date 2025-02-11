**SSRseq_count: A Next-Generation Sequencing (NGS)-Based Microsatellite (SSR) Genotyping Tool**  
*An optimized version of [ccoo22/SSRseq_count](https://github.com/ccoo22/SSRseq_count) for enhanced accuracy and compatibility.*  

---

### Overview  
This repository provides a modified version of the SSRseq_count tool for microsatellite (SSR) genotyping using high-throughput sequencing data. The original tool, designed for polyploid species, resolves allele dosage uncertainty and improves analyses of genetic diversity, population structure, and differentiation (see [References](#references)). This fork introduces critical optimizations for paired-end read processing and motif determination, ensuring compatibility with modern bioinformatics workflows.  

---

### Key Modifications  
1. **Paired-End Read Merging**  
   - **Overlapping reads**: Merged using [FLASH v1.2.11](https://ccb.jhu.edu/software/FLASH/) for high-accuracy overlap assembly.  
   - **Non-overlapping reads**: Directly concatenated to retain maximum sequence information.  

2. **Genotyping Workflow**  
   - Replaced Perl-based regular expression motif matching with **BLAST alignment** for robust and accurate motif determination.  
   - Ensures compatibility with updated BLAST versions (fixes errors in legacy implementations).  

---

### Script Versions  
| File | Description |  
|------|-------------|  
| [`str_count_orig.pl`](str_count_orig.pl) | Original script from [ccoo22/SSRseq_count](https://github.com/ccoo22/SSRseq_count). |  
| [`str_count_reg.pl`](str_count_reg.pl) | Updated version with BLAST command fixes for compatibility with modern BLAST versions. |  
| [`str_count_ovlp.pl`](str_count_ovlp.pl) | Processes **only overlapping paired-end reads** (uses FLASH-merged reads). |  
| [`str_count.pl`](str_count.pl) | Processes **overlapping paired-end reads** (FLASH-merged) and **non-overlapping reads**. Implements **BLAST alignment** (replacing Perl regex motif matching) for enhanced accuracy in motif determination. |

---

### References  
Cui X, Li C, Qin S, et al. High-throughput sequencing-based microsatellite genotyping for polyploids to resolve allele dosage uncertainty and improve analyses of genetic diversity, structure and differentiation: A case study of the hexaploid *Camellia oleifera*. *Mol Ecol Resour*. 2021;00:1–14. [doi:10.1111/1755-0998.13469](https://doi.org/10.1111/1755-0998.13469). PMID: [34260828](https://pubmed.ncbi.nlm.nih.gov/34260828/).  

---

### Usage Notes  
- **Dependencies**: Ensure FLASH v1.2.11 and BLAST are installed.  
- **Input**: Paired-end FASTQ files and reference SSR motifs.  
- **Output**: Allele counts and genotypes compatible with downstream population genetics analyses.  

For detailed instructions, refer to the original repository or contact contributors.  

---  
🔬 Built for polyploid genomics | 🛠 Optimized for accuracy and scalability | 💡 Open-source under MIT License

