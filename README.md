[README (2).md](https://github.com/user-attachments/files/25162750/README.2.md)
# BCR Analysis Suite - Integrated CompIgS + FastBCR Analysis Platform

<img width="265" height="112" alt="BCR Analysis Suite Logo" src="https://github.com/user-attachments/assets/c7e5c63b-148f-4f36-b121-344da1aeeea9" />

**BCR Analysis Suite** is an open-source computational platform that combines two powerful methods for comparative analysis of B cell receptor (BCR) repertoires across antibody isotypes.

---

## 📅 Introduction

Antibodies, or immunoglobulins (Ig), are produced by B cells and form a diverse repertoire critical for adaptive immunity. Each antibody is defined by unique V(D)J gene recombination forming clonotypes, which may exist in various isotypes (IgM, IgE, IgG, IgA, IgD, etc.) with different effector functions.

Affinity maturation and somatic hypermutation refine antibody specificity over time. Understanding how clonotypes are shared or diverge between different antibody isotypes provides critical insights into immune responses.

**BCR Analysis Suite** integrates two complementary approaches:

1. **CompIgS Method** - Uses exact V+D+J gene matching to identify clonotypes
2. **FastBCR Method** - Uses V+J genes with k-mer seeding and junction similarity analysis

This dual-method approach enables comprehensive analysis of clonal relationships, somatic hypermutation patterns, and divergent clone identification across antibody isotypes.

---

## 🔬 Analysis Methods

### CompIgS Method (Python/VDJ Matching)

**Input Requirements:**
- IMGT/HighV-QUEST output files
- IMGT/StatClonotype files
- File types: `.txt` files containing:
  - StatClonotype data (both isotypes)
  - V-REGION-nt-mutation-statistics (both isotypes)
  - AA-sequences (both isotypes)
  - V-REGION-mutation-and-AA-change-table (both isotypes)
  - IMGT VH reference sequences (FASTA)

**Method:**
- Identifies clonotypes based on **exact V+D+J gene matching**
- Analyzes shared and unique clones between isotypes
- Calculates somatic hypermutation (SHM) rates
- Identifies divergent clones based on CDR3 amino acid differences
- Categorizes mutations as low (0-5), moderate (6-11), or high (≥12)

**Key Features:**
- Clonal expansion analysis
- Copy number ratio calculations between isotypes
- Mutation profile comparison (CDR2 and CDR3 regions)
- VH sequence extraction and completion
- Amino acid composition analysis using EMBOSS categories

### FastBCR Method (R/VJ + K-mer + Junction Similarity)

**Input Requirements:**
- AIRR-formatted TSV files (`db-pass.tsv` from Change-O/IgBLAST pipeline)
- Files needed:
  - Isotype 1 TSV file (e.g., IgE)
  - Isotype 2 TSV file (e.g., IgG1)

**Method:**
- **Within-isotype clustering:** Uses V+J genes and k-mer seeding (six 5-mers) with consensus scoring (≥0.8)
- **Between-isotype sharing:** Identifies shared clusters using:
  - V+J gene matching
  - Configurable junction amino acid similarity threshold (default 80%)
- Minimum amino acid difference threshold for divergent clones (default ≥2)
- Minimum cluster depth filtering (default ≥3 sequences)

**Key Features:**
- Hierarchical clustering within each isotype
- Shared cluster identification between isotypes
- Divergent clone detection at all mutation levels
- SHM analysis for shared and divergent populations
- Lineage tree visualization for clonal relationships
- Multiple sequence alignments (MSA) for divergent sequences
- Comprehensive statistical reporting

**Key Difference from CompIgS:**
- FastBCR does not require D-gene matching and uses junction sequence similarity
- More flexible for datasets where D-gene assignment may be ambiguous
- Optimized for parallel processing (Windows and Unix compatible)

---

## 📦 Key Features

* 🔍 **Dual Analysis Methods**: Choose between exact VDJ matching (CompIgS) or VJ + junction similarity (FastBCR)
* 🧬 **Comprehensive BCR Profiling**: Analyze clonal expansion, mutation patterns, and divergence
* 📊 **Rich Visualizations**: Automated generation of plots, heatmaps, lineage trees, and MSA alignments
* 🧠 **Mutation Level Analysis**: Categorize clones as low, moderate, or highly mutated
* 🔬 **Divergent Clone Detection**: Identify isotype-specific divergent sequences with configurable thresholds
* 🖥️ **GUI Interface**: User-friendly graphical interface for both methods (no coding required)
* 🪟 **Cross-Platform**: Windows (.exe), macOS (.app), and Python/R versions available
* 🚀 **Performance Optimized**: Parallel processing support for large datasets
* 📈 **Export Options**: Multiple output formats including CSV, PDF, PNG, and FASTA

---

## 🔹 Table of Contents

* [Quick Start](#-quick-start)
* [Installation](#-installation)
* [Input Data Requirements](#-input-data-requirements)
* [Usage Examples](#-usage-examples)
* [Output Files](#-output-files)
* [Analysis Parameters](#-analysis-parameters)
* [Citation](#-citation)
* [License](#-license)
* [Contact](#-contact)

---

## 🚀 Quick Start

### Option 1: Windows Executable (No Installation Required)

Download the standalone `.exe` from:
* 🔗 [GitHub Releases](https://github.com/Chrisjames1992/CompIgS/releases)
* 📓 [Zenodo DOI](https://doi.org/10.5281/zenodo.15774119)

Double-click to run - no Python or R installation needed!

### Option 2: macOS Application

Download the macOS version:
* 🔗 [GitHub Releases](https://github.com/Chrisjames1992/CompIgS/releases/download/mv.01/CompIgS_GUI.app.zip)

### Option 3: Python GUI

```bash
# Clone the repository
git clone https://github.com/Chrisjames1992/CompIgS.git
cd CompIgS

# Run the GUI
python CompIgS_FastBCR.py
```

---

## 🛠️ Installation

### Pre-built Applications (Recommended)

**Windows Executable (.exe):**
- Download from [GitHub Releases](https://github.com/Chrisjames1992/CompIgS/releases)
- Download from [Zenodo](https://doi.org/10.5281/zenodo.15774119)
- No installation required - just download and run!

**macOS Application (.app):**
- Download from [GitHub Releases](https://github.com/Chrisjames1992/CompIgS/releases/download/mv.01/CompIgS_GUI.app.zip)
- Unzip and drag to Applications folder

**Python GUI:**
- Clone repository: `git clone https://github.com/Chrisjames1992/CompIgS.git`
- Run: `python CompIgS_FastBCR.py`
- All Python and R code included in the repository

### System Requirements

- **RAM:** 8GB minimum, 16GB+ recommended for large datasets
- **Storage:** 1GB+ free space for output files
- **OS:** Windows 10+, macOS 10.14+, or Linux (Ubuntu 18.04+)
- **Python:** ≥ 3.8 (for Python GUI)
- **R:** ≥ 4.0 (for FastBCR method)

---

## 📂 Input Data Requirements

### CompIgS Method Input Files

Import the following **9 IMGT files** (use Isotype 1 as your query isotype, e.g., IgE):

1. `stats_IgE.txt` - IgE IMGT StatClonotype
2. `stats_IgG1.txt` - IgG1 IMGT StatClonotype
3. `stats_IgE_mut.txt` - IgE\_8\_V-REGION-nt-mutation-statistics
4. `stats_IgG1_mut.txt` - IgG1\_8\_V-REGION-nt-mutation-statistics
5. `stats_IgE_Aminotab.txt` - IgE\_5\_AA-sequences
6. `stats_IgG1_Aminotab.txt` - IgG1\_5\_AA-sequences
7. `stats_IgE_Aminotab_change.txt` - IgE\_7\_V-REGION-mutation-and-AA-change-table
8. `stats_IgG1_Aminotab_change.txt` - IgG1\_7\_V-REGION-mutation-and-AA-change-table
9. `sequences.fasta` - Human or Mouse VH reference (IMGT format)

**Required IMGT Columns:**
- `V_gene`, `D_gene`, `J_gene`
- `Sequence.ID` or `seqid`
- `V-DOMAIN Functionality`
- `V-REGION Nb of nucleotides`
- `CDR3-IMGT Nb of nonsilent mutations`
- `onecopy` (copy number/clone size)

### FastBCR Method Input Files

Requires **2 AIRR-formatted TSV files** (output from Change-O or IgBLAST):

1. Isotype 1 file (e.g., `IgE_db-pass.tsv`)
2. Isotype 2 file (e.g., `IgG1_db-pass.tsv`)

**Required AIRR Columns:**
- `v_call`, `j_call` (V and J gene assignments)
- `junction_aa` (CDR3 junction amino acid sequence)
- `sequence_id` (unique sequence identifier)
- `productive` (TRUE/FALSE for productive rearrangements)
- Optional: `sequence_alignment`, `germline_alignment` (for SHM calculation)

**Accepted File Formats:**
- `.tsv` (tab-separated values)
- Must be AIRR-compliant format from Change-O pipeline

---

## 💻 Usage Examples

### GUI Usage (Recommended)

1. **Launch the application**
   - Windows: Double-click `CompIgS_FastBCR.exe`
   - macOS: Open `CompIgS_GUI.app`
   - Python: Run `python CompIgS_FastBCR.py`

2. **Select Analysis Method**
   - Choose "CompIgS" tab for VDJ exact matching
   - Choose "FastBCR" tab for VJ + junction similarity analysis

3. **Configure Analysis**
   - Set isotype names (e.g., IgE, IgG1)
   - Browse and select input files
   - Adjust analysis parameters (thresholds, filters)
   - Choose output directory

4. **Run Analysis**
   - Click "Start Analysis"
   - Monitor progress in the log viewer
   - Review results when complete

### Command Line Usage (Python)

```python
from BCR_Analysis_Suite import CompIgSAnalyzer, FastBCRAnalyzer

# CompIgS Analysis
analyzer = CompIgSAnalyzer(
    iso1_name="IgE",
    iso2_name="IgG1",
    cdr3_threshold=2,
    cdr3cdr2_threshold=2
)

analyzer.load_files(
    stats_iso1="data/stats_IgE.txt",
    stats_iso2="data/stats_IgG1.txt",
    # ... additional files
)

results = analyzer.run_analysis()

# FastBCR Analysis
fastbcr = FastBCRAnalyzer(
    iso1_file="data/IgE_db-pass.tsv",
    iso2_file="data/IgG1_db-pass.tsv",
    output_dir="results/",
    junction_threshold=80,
    min_aa_diff=2,
    min_depth=3
)

fastbcr.run()
```

### R Script Usage (FastBCR)

```r
# Load required libraries
library(fastBCR)
library(ggplot2)
library(dplyr)

# Set parameters
ISO1_FILE <- "data/IgE_db-pass.tsv"
ISO2_FILE <- "data/IgG1_db-pass.tsv"
JUNCTION_SIMILARITY_THRESHOLD <- 80
DIVERGENT_MIN_DIFFERENCES <- 2

# Run analysis
source("fastBCR_analysis_script.R")
```

---

## 📄 Output Files

### CompIgS Outputs

**Generated in `~/Downloads/CompIgS_Analysis/`:**

**CSV Files (18 total):**
- `captured.csv` - Summary of all captured clones
- Isotype-specific population files (17 files):
  - Total, shared, unique populations
  - Mutated/unmutated subsets
  - Biased clones (Isotype 1-biased, Isotype 2-biased)
  - Divergent clones by mutation level

**Plots Folder (`BCR_Analysis_Plots/`):**
- Clone distribution plots
- Copy number ratio visualizations
- Mutation distribution graphs
- Clonal size distribution plots (zoomed views)
- CDR2/CDR3 mutation comparison bar plots
- Mutated vs unmutated pie charts

**VH Sequences Folder:**
- Completed VH sequences for divergent clones
- Top 10 expanded clones (by mutation level)
- Clonally related sequences

**Data Folder:**
- Amino acid composition analysis (EMBOSS categories)
- Excel files for Isotype 1 and Isotype 2

### FastBCR Outputs

**Generated in specified output directory:**

**CSV Files:**
- `IgE_cluster_SHM.csv` - SHM statistics per cluster (Isotype 1)
- `IgG1_cluster_SHM.csv` - SHM statistics per cluster (Isotype 2)
- `IgE_IgG1_shared_clusters.csv` - Shared cluster details
- `Shared_IgE_clusters_SHM.csv` - SHM for shared Isotype 1 clusters
- `Divergent_IgE_clones_all_levels.csv` - All divergent clones
- `Divergent_IgE_clones_filtered.csv` - Filtered divergent clones (≥2 per cluster)
- `NonDivergent_IgE_by_mutation_level.csv` - Non-divergent statistics
- `IgE_IgG1_summary_statistics.csv` - Complete analysis summary

**Visualization Files (PDF + PNG):**

*Cluster Analysis:*
- `IgE_clusters_bubble.pdf/png` - Bubble plot of clusters
- `IgE_V_gene_usage.pdf/png` - V gene usage pie chart
- `IgE_J_gene_usage.pdf/png` - J gene usage pie chart
- `IgE_VJ_pairing.pdf/png` - V-J pairing heatmap
- `IgE_junction_length.pdf/png` - Junction length distribution

*Shared Cluster Analysis:*
- `Shared_clusters_multi_panel.pdf/png` - 4-panel shared cluster overview
- `Shared_clusters_heatmap.pdf/png` - Cluster matching heatmap
- `Junction_similarity_distribution.pdf/png` - Junction similarity histogram
- `VJ_combinations_shared_clusters.pdf/png` - Top V/J combinations
- `Cluster_overlap_pie_charts.pdf/png` - Overlap percentages
- `Shared_IgE_SHM_ascending.pdf/png` - SHM per cluster

*Divergent Clone Analysis:*
- `Divergent_IgE_by_mutation_level.pdf/png` - Clusters with divergent clones
- `Total_divergent_IgE_by_mutation_level.pdf/png` - Total divergent sequences
- `Divergent_IgE_SHM_vs_pct.pdf/png` - SHM correlation
- `Divergent_IgE_copy_sizes.pdf/png` - Copy size comparison
- `Divergent_IgE_vs_IgG1_sizes.pdf/png` - Scatter plot comparison
- `Divergent_IgE_summary_table.pdf/png` - Summary table

*Lineage Trees (Top N clusters):*
- `Shared_RankN_IgE{id}_IgG1{id}_lineage_tree.pdf/png`
- `Divergent_{Level}_RankN_IgE{id}_IgG1{id}_lineage_tree.pdf/png`
- Corresponding `.nwk` files (Newick format)

*Multiple Sequence Alignments:*
- `Shared_RankN_IgE{id}_IgG1{id}_alignment.fasta`
- `Shared_RankN_IgE{id}_IgG1{id}_MSA.pdf/png`
- `Shared_RankN_IgE{id}_IgG1{id}_MSA_pretty.pdf` (formatted MSA)

**Log Files:**
- `session_info.txt` - R session information
- Performance statistics and memory usage

---

## ⚙️ Analysis Parameters

### CompIgS Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cdr3_threshold` | 2 | Minimum CDR3 AA differences for divergent clones |
| `cdr3cdr2_threshold` | 2 | Minimum CDR3+CDR2 AA differences for divergence |
| `productive_only` | True | Filter for productive sequences only |
| `min_copy_number` | 2 | Minimum clone copy number for inclusion |

**Mutation Level Thresholds:**
- Low: 0-5 mutations
- Moderate: 6-11 mutations
- High: ≥12 mutations

### FastBCR Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `junction_threshold` | 80 | Junction AA similarity threshold (%) for shared clusters |
| `min_aa_diff` | 2 | Minimum AA differences for divergent classification |
| `min_depth` | 3 | Minimum sequences per cluster |
| `max_depth` | 1000 | Maximum cluster size threshold |
| `overlap_threshold` | 0.1 | K-mer overlap threshold for clustering |
| `consensus_threshold` | 0.80 | Consensus score for cluster assignment |
| `productive_only` | False | Filter productive sequences (optional) |
| `n_top_shared` | 10 | Number of top shared clusters for trees/MSA |
| `n_top_divergent` | 10 | Number of top divergent clusters per mutation level |

**SHM Calculation Methods (Auto-detected):**
- Alignment-based (IMGT style)
- Pre-calculated mutation rate (MiXCR style)
- V-identity based (derived from v_identity field)

---

## 🧪 Example Workflow

### Complete Analysis Pipeline

1. **Data Preparation**
   ```bash
   # Process sequences through IMGT/HighV-QUEST or Change-O
   # Organize output files into input directory
   ```

2. **CompIgS Analysis (VDJ Exact Matching)**
   ```bash
   # Launch GUI and select CompIgS tab
   # Load all 9 IMGT files
   # Configure parameters
   # Run analysis
   # Review 18 CSV outputs and visualizations
   ```

3. **FastBCR Analysis (VJ + Junction Similarity)**
   ```bash
   # Launch GUI and select FastBCR tab
   # Load 2 AIRR TSV files
   # Set junction similarity threshold (80%)
   # Configure divergent clone parameters
   # Run analysis
   # Review clusters, SHM, and lineage trees
   ```

4. **Compare Results**
   - CompIgS: Strict VDJ matching, precise clonotype definition
   - FastBCR: Flexible VJ matching, junction-based similarity
   - Identify consensus clones found by both methods
   - Investigate method-specific clones

---

## 📊 Interpretation Guide

### Understanding Shared vs Divergent Clones

**Shared Clones:**
- Found in both isotypes with same V(D)J genes (CompIgS) or V+J genes (FastBCR)
- Represent class-switched B cells
- May show different mutation levels between isotypes

**Divergent Clones (Isotype 1-specific):**
- Present in Isotype 1 (e.g., IgE) shared cluster
- Junction sequence differs by ≥N amino acids from all Isotype 2 sequences
- Indicates isotype-specific somatic hypermutation
- Critical for understanding allergen-specific responses

**Mutation Level Categories:**
- **Low (0-5 mutations):** Recently activated or germline-like
- **Moderate (6-11 mutations):** Undergoing affinity maturation
- **High (≥12 mutations):** Highly mutated, affinity-matured clones

### Key Metrics

**Copy Number Ratio (Isotype 1 / Isotype 2):**
- Ratio > 2: Isotype 1-biased expansion
- Ratio 0.5-2: Balanced expansion
- Ratio < 0.5: Isotype 2-biased expansion

**Junction Similarity (FastBCR):**
- ≥80%: Likely clonally related
- 60-79%: Possibly related (review V/J genes)
- <60%: Unlikely to be clonally related

---

## 📚 Citation

If you use **BCR Analysis Suite** or **CompIgS**, please cite:

> Udoye, C. C., Mehrabani Khasraghi, S., Witt, P., Biswas, S., Manz, R., & Fähnrich, A. (2025). *CompIgS: A Computational Workflow for Comparative Analysis of Related Clonotypes within Distinct Antibody Subclasses*. BMC Bioinformatics.

If you use **FastBCR** specifically, please also cite:

> Wang, K., & Zhang, J. (2023). *fastBCR: a k-mer-based approach for fast inferencing clonal families from large-scale B cell repertoire sequencing data*. Briefings in Bioinformatics, 24(6), bbad404. https://doi.org/10.1093/bib/bbad404
> 
> **PMC Article:** https://pmc.ncbi.nlm.nih.gov/articles/PMC10626204/
> 
> **GitHub Repository:** https://github.com/ZhangLabTJU/fastBCR

---

## 📖 Resources

### Publications

**CompIgS:**
- Udoye, C. C., et al. (2025). BMC Bioinformatics.

**FastBCR:**
- Wang, K., & Zhang, J. (2023). *fastBCR: a k-mer-based approach for fast inferencing clonal families from large-scale B cell repertoire sequencing data*. Briefings in Bioinformatics, 24(6), bbad404.
- DOI: https://doi.org/10.1093/bib/bbad404
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10626204/

### GitHub Repositories

**BCR Analysis Suite (CompIgS + FastBCR):**
- https://github.com/Chrisjames1992/CompIgS

**Original FastBCR:**
- https://github.com/ZhangLabTJU/fastBCR

**IMGT/HighV-QUEST (for data generation):**
- http://www.imgt.org/HighV-QUEST/

---

## 🖼️ Logo

You may use the following logo with attribution in presentations and publications:

<img width="265" height="112" alt="BCR Analysis Suite Logo" src="https://github.com/user-attachments/assets/c7e5c63b-148f-4f36-b121-344da1aeeea9" />

---

## 📜 License

MIT License – see `LICENSE` for details.

---

## 🙋‍♀️ Contact & Support

For questions, suggestions, or bug reports:

* 📧 Open an issue on [GitHub](https://github.com/Chrisjames1992/CompIgS/issues)
* 💬 Discussions: Use GitHub Discussions for methodology questions
* 📖 Documentation: See [Wiki](https://github.com/Chrisjames1992/CompIgS/wiki)

---

## 🔄 Version History

**v2.0.0** (Current)
- Integrated FastBCR method alongside CompIgS
- Added dual-method GUI interface
- Enhanced divergent clone detection
- Lineage tree and MSA visualization
- Parallel processing optimization

**v1.0.0**
- Initial CompIgS release
- VDJ exact matching
- Basic mutation analysis

---

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

Areas where contributions are especially valuable:
- Additional visualization options
- Support for other BCR analysis tools
- Performance optimizations
- Documentation improvements
- Bug reports and fixes

---

## ⚠️ Important Notes

1. **Data Privacy**: All analysis is performed locally - no data is sent to external servers
2. **Memory Requirements**: Large datasets (>100K sequences) may require 16GB+ RAM
3. **Processing Time**: FastBCR analysis with lineage trees can take 30-60 minutes for large datasets
4. **File Size**: IMGT output files can be large; ensure adequate storage space
5. **Column Names**: The tool auto-detects column name variations (e.g., "Sequence ID" vs "Sequence.ID")

---

## 🔍 Troubleshooting

**Common Issues:**

1. **"R not found" error (FastBCR)**
   - Install R from https://cran.r-project.org/
   - Ensure R is in system PATH
   - Restart application after R installation

2. **"fastBCR package missing"**
   - The fastBCR R package and all dependencies are included in the repository
   - See https://github.com/ZhangLabTJU/fastBCR for package documentation

3. **Out of memory errors**
   - Reduce dataset size or increase system RAM
   - Close other applications
   - Use productive-only filtering

4. **Empty results**
   - Check input file format (IMGT vs AIRR)
   - Verify column names match expected format
   - Ensure productive sequences present
   - Review filtering thresholds

5. **Visualization errors**
   - Ensure all required R packages installed
   - Check output directory write permissions
   - Review log for specific error messages

For additional help, see [Troubleshooting Guide](https://github.com/Chrisjames1992/CompIgS/wiki/Troubleshooting).

---

**Built with ❤️ for the BCR analysis community**
