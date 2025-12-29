# Compare-Gene-Expressions
A Python-based computational framework for analyzing differential gene expression using GEO microarray data.
## 🚀 Features
- **Smart Metadata Parser:** Automatically identifies tissue groups (Brain, Lung, Heart, etc.) from Series Matrix files.
- **Interactive Search:** Find specific Gene/Probe IDs without opening large text files.
- **Statistical Analysis:** Performs Welch’s T-test to determine differential expression.
- **Visualization:** Generates publication-ready boxplots using Seaborn.

## 🛠️ Installation & Usage
1. Clone the repo: `git clone https://github.com/YOUR_USERNAME/BioAnalyzer-GEO.git`
2. Install dependencies: `pip install pandas seaborn scipy matplotlib`
3. Instiall the GEO microarray data from `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/matrix/`
4. Place `GSE1000_series_matrix.txt` in the root folder.
5. Run the tool: `python analysis.py`

## 📊 Example Output
The tool identifies tissue differences with high precision. For example, comparing GAPDH (201210_at) across organs provides a baseline for housekeeping gene stability.
