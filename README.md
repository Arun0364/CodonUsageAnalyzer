# 🧬 Codon Usage Analyzer

A Python-based bioinformatics tool that analyzes **codon usage frequencies** in DNA coding sequences (CDS).  
The script counts codons, calculates usage bias, and generates **bar plots, heatmaps, and comparison plots** between species or genes.  

---

## ✨ Features
- Accepts **FASTA DNA sequences** as input  
- Counts codon frequencies and converts them to relative usage (%)  
- Generates:
  - Codon usage barplots (per sequence)
  - Codon usage heatmaps (codon table layout)
  - Side-by-side barplots for two sequences
  - Difference heatmaps showing codon bias shifts
- Exports results to CSV and PNG  

---

## 📂 Repository Structure
```CodonUsageAnalyzer/
├── codon_usage.py          # Main Codon Usage script 
├── README.md               # Documentation for the project
├── outputs                 # Sample output heatmap of the alignment matrix               
└── sequences/              # Folder containing sample input sequences
    ├── Human_HBB.fasta
    └── Mouse_HBB.fasta
├── LICENSE
└── .gitignore
```

---

## 🚀 Usage

### 1. Install dependencies
```bash
pip install -r requirements.txt
```
### 2. Run the program
```bash
python3 codon_usage.py
```

### 3. Provide input FASTA files
You will be prompted for:
First FASTA file path + label (e.g., Human_HBB.fasta, "Human")
Second FASTA file path + label (optional, e.g., Mouse_HBB.fasta, "Mouse")

---

## 📊 Example Outputs

+ codon_usage_Human.csv → raw counts + frequency table
+ codon_plot_Human.png → codon barplot
+ codon_heatmap_Human.png → codon heatmap
+ codon_comparison.png → side-by-side barplot (if 2 sequences given)
+ codon_diff_heatmap.png → codon bias difference heatmap

---

## 🛠️ Requirements
* Python 3.8+
* Biopython
* Pandas
* Matplotlib
* Seaborn

Install via:

```bash
pip3 install biopython pandas seaborn matplotlib

```

---

## ⚙️ How It Works

The workflow for codon usage analysis:
``` mermaid
flowchart TD
    A[📥 Input FASTA file(s)] --> B[🔍 Parse DNA sequence]
    B --> C[🧮 Count codons (triplets)]
    C --> D[📊 Calculate codon frequencies (%)]
    D --> E[💾 Export CSV results]
    D --> F[📈 Visualizations]
    F --> F1[Barplots per sequence]
    F --> F2[Heatmaps per sequence]
    F --> F3[Comparison plots (2 sequences)]
    F --> F4[Difference heatmap (bias shift)]
    E --> G[📂 Outputs directory]
    F1 --> G
    F2 --> G
    F3 --> G
    F4 --> G
```

---

## 🏷️ Tags

`bioinformatics` `genomics` `codon-usage` `sequence-analysis` `computational-biology` `python` `data-visualization` `heatmap` `biopython` `molecular-biology` `FASTA`

---

## 🤝 Contributing

Pull requests are welcome. If you’d like to extend this project (e.g., codon optimization, RSCU index, or multi-species batch analysis), feel free to fork and contribute.

---

## 👤 Contact

Arunannamalai Sujatha Bharath Raj

📧 [arun03bt@gmail.com]

🔗 [https://www.linkedin.com/in/arunannamalai-sb-823351344/](https://www.linkedin.com/in/arun-823351344/)

🐙 [https://github.com/Arun0364](https://github.com/Arun0364)

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## ❤️ Acknowledgements

Biopython – for FASTA parsing and sequence utilities

Pandas – for tabular data handling and CSV exports

Seaborn & Matplotlib – for barplots and heatmaps

Genomic Data Repositories – NCBI and Ensembl for CDS FASTA files

Codon Bias Research – decades of work on codon usage in molecular biology and comparative genomics
