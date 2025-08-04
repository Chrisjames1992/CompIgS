# CompIgS – Comparative Analysis of Clonotypes Across Antibody Subclasses

<img width="265" height="112" alt="image" src="https://github.com/user-attachments/assets/c7e5c63b-148f-4f36-b121-344da1aeeea9" />


**CompIgS** is an open-source computational tool for the comparative analysis of related clonotypes across antibody subclasses from the same sample.

---

## 📅 Introduction

Antibodies, or immunoglobulins (Ig), are produced by B cells and form a diverse repertoire critical for adaptive immunity. Each antibody is defined by unique V(D)J gene recombination forming clonotypes, which may exist in various subclasses. These subclasses, such as IgM, IgE, IgG, etc., have different effector functions.

Affinity maturation and somatic hypermutation refine antibody specificity. However, tools for comparing clonotypes across isotypes remain limited.

**CompIgS** fills this gap by identifying shared and divergent clonotypes using annotated output from IMGT/HighV-QUEST and IMGT/StatClonotype. It facilitates clonal expansion, mutation profile comparison, and divergence analysis between isotypes such as IgE and IgG1, as demonstrated in a murine food allergy model.

---

## 📦 Key Features

* 🔍 Compare antibody isotypes (e.g., IgE (Ig1) vs. IgG1 (Ig2)) in one sample
* 🧠 Filter clonotypes based on various metrics: divergent CDR3AA, clonal size, shared vs. unique, mutated vs. unmutated
* 🖥️ GUI for local, interactive analysis
* ☁️ Google Colab notebook for quick, cloud-based access
* 🫿 Windows executable – no coding required
* Mac Os version-no coding required

---

## 🔹 Table of Contents

* [Quick Start](#-quick-start)
* [Overview](#-project-structure)
* [Data/Samples](#-input-data-requirements)

---

## 🚀 Quick Start

### Get the Software

**Option 1: Run the Windows `.exe` version**
Download the standalone executable from:

* 🔗 [GitHub Releases](https://github.com/Chrisjames1992/CompIgS/releases)
* 📓 [Zenodo DOI](https://doi.org/10.5281/zenodo.15774119)

**Option 2: Run the macOS version**
Download the standalone from:
* 🔗 [GitHub Releases](https://github.com/Chrisjames1992/CompIgS/releases)

**Option 3: Use via Google Colab**
Open and run `CompIgS.ipynb` in Google Colab and follow the instructions to upload your data and begin analysis.

**Option 4: Run GUI locally (Python required)**

### 🛠️ Requirements

* Python ≥ 3.1
* Libraries: `pandas`, `numpy`, `PyQt5`, `seaborn`, `matplotlib`, `peptides`

---

## 📁 Project Structure

```
CompIgS/
├── CompIgS.ipynb         # Google Colab notebook
├── CompIgS_GUI.py        # GUI version for local use
├── README.md             # This file
└── releases/             # Downloadable Windows and macOS EXE binaries
```

---

## 📂 CompIgS Notebook: Local GUI and Colab Versions

We offer one notebook version and one `.py` file, tailored for cloud-based (Colab) execution or local (GUI).

### Google Colab Version

* `CompIgS.ipynb` (without GUI)
* Usage: Upload your IMGT/HighV-QUEST and IMGT/StatClonotype output data and follow the instructions
* Open notebook in Google Colab
* Dependencies: `pandas`, `numpy`, `biopython`, `seaborn`, `matplotlib`, `google.colab`, `collections`
* Example data: available on Google Colab or Zenodo

### Input Data Requirements (for both versions)

Import the following IMGT/StatClonotype and IMGT/HighV-QUEST .txt files in the directory or use example_data [example data](https://doi.org/10.5281/zenodo.15774119) . Use Ig1 for the query isotype, in the example IgE is Ig1:

* `stats_IgE.txt` - IgE IMGT StatClonotype
* `stats_IgG1.txt` - IgG1 IMGT StatClonotype
* `stats_IgE_mut.txt` - IgE\_8\_V-REGION-nt-mutation-statistics
* `stats_IgG1_mut.txt` - IgG1\_8\_V-REGION-nt-mutation-statistics
* `stats_IgE_Aminotab.txt` - IgE\_5\_AA-sequences
* `stats_IgG1_Aminotab.txt` - IgG1\_5\_AA-sequences
* `stats_IgE_Aminotab_change.txt` - IgE\_7\_V-REGION-mutation-and-AA-change-table
* `stats_IgG1_Aminotab_change.txt` - IgG1\_7\_V-REGION-mutation-and-AA-change-table
* `sequences.fasta` - Human or Mouse VH reference

---

## 🧪 Input Format

* Accepted formats: `.txt` files from IMGT/StatClonotype and IMGT/HighV-QUEST
* Required fields: `V_gene`, `D_gene`, `J_gene`, `Sequence.ID`, `V-DOMAIN Functionality`, `V-REGION Nb of nucleotides`, `CDR3-IMGT Nb of nonsilent mutations`, etc.

---

## 📄 Output

* Summary plots saved in a folder named `BCR_Analysis_Plots` in your Downloads directory
* 18 CSV files: one `captured.csv` and 17 describing different Ig1 and Ig2 subpopulations

---

## 📚 Citation

If you use **CompIgS**, please cite:

> Udoye, C. C., Mehrabani Khasraghi, S., Witt, P., Biswas, S., Manz, R., & Fähnrich, A. (2025). *CompIgS: A Computational Workflow for Comparative Analysis of Related Clonotypes within Distinct Antibody Subclasses*. BMC Bioinformatics.

---

## 🖼️ Logo

You may use the following logo with attribution in presentations and publications:

<img width="265" height="112" alt="image" src="https://github.com/user-attachments/assets/c7e5c63b-148f-4f36-b121-344da1aeeea9" />

---

## 📜 License

MIT License – see `LICENSE` for details.

---

## 🙋‍♀️ Contact

For questions, suggestions, or bug reports:

* Open an issue on [GitHub](https://github.com/Chrisjames1992/CompIgS/issues)
