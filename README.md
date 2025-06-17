# Genetic landscape of an _in vivo_ protein interactome

## __Université de Montréal (2021-2025)__

- Savandara Besse __#__, Tatsuya Sakaguchi __#__, Louis Gauthier, Zahra Sahaf, Olivier Péloquin, Lidice Gonzalez, Xavier Castellanos-Girouard, Nazli Koçatug, Chloé Matta, Julie Hussin, Stephen Michnick\*, Adrian Serohijos\*
  - Department of Biochemistry, Université de Montréal, Montréal, Québec, Canada (SB, TS, LoG, ZS, OP, LiG, XCG, NK, CM, SWM, AS)
  - Robert-Cedergren Center for Bioinformatics and Genomics, Université de Montréal, Montréal, Québec, Canada (SB, TS, LoG, ZS, OP, LiG, XCG, NK, CM, SM, AS)
  - Institut de Cardiologie de Montréal, Montréal, Québec, Canada (JH)
  - Département de Médecine, Faculté de Médecine, Université de Montréal, Montréal, Québec, Canada (JH)

\# These authors contributed equally.

## piQTL website

> https://ladyson1806.github.io/SerohijosLab-piQTL/

## 1. Contact

### \* Correspondence

- Adrian Serohijos (@aserohijos) : adrian.serohijos @ umontreal.ca
- Stephen Michnick (@michnics) : stephen.michnick @ umontreal.ca

### Website and Code maintenance

- Savandara Besse (@ladyson1806) : savandara.besse @ cnrs.fr


## 2. Software and libraries

### Python

- Python 3.11.4
- [matplotlib 3.7.2](https://pypi.org/project/matplotlib/)
- [matplotlib-venn 0.11.9](https://pypi.org/project/matplotlib-venn/)
- [mpl-scatter-density 0.7](https://pypi.org/project/mpl-scatter-density/)
- [numpy 1.25.2](https://pypi.org/project/numpy/)
- [pandas 2.0.3](https://pypi.org/project/pandas/)
- [plotnine 0.12.2](https://pypi.org/project/plotnine/)
- [scipy 1.11.1](https://pypi.org/project/scipy/)
- [seaborn 0.12.2](https://pypi.org/project/seaborn/)
- [tqdm 4.66.1](https://pypi.org/project/tqdm/)
- [venn 0.1.3](https://pypi.org/project/venn/)

### R

- R (ToDo: Specify the R version in this study)
- Required R packages: (ToDo: please specify as needed)

---

## 3. Installation Guide

### Python

1. **Clone the repository:**
    ```bash
    git clone https://github.com/ladyson1806/SerohijosLab-piQTL.git
    cd SerohijosLab-piQTL
    ```
    *Note: The repository size is approximately 15 GB. Depending on your internet speed and disk performance, cloning may take 10–30 minutes or longer.*

2. **(Recommended) Create a virtual environment:**
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```

3. **Install Python dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

### R

1. **Install R (version 4.2 or higher recommended).**

2. **Install required R packages:**
    ```R
    install.packages(c("tidyverse", "data.table", ...)) # Add all required packages here
    ```

---

## 4. Instructions for Use
The sequence data generated in this study are available in the Gene Expression Omnibus (GEO) under accession ID [GSE246414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246414).
All Python scripts, R scripts, and Jupyter notebooks used for data analysis are provided and are fully executable.
The outputs from both the analysis scripts and notebooks, including all figures, are saved in the `data/` and `figure/` directories.
To reproduce these results, execute the scripts and notebooks as described belolw.

### Running Python Scripts

To run a Python analysis script:
```bash
python target_script.py
```
Replace `target_script.py` with the name of the script you wish to execute.

### Running R Scripts

To run an R script:
```bash
Rscript target_script.R
```
Replace `target_script.R` with the name of the script you wish to execute.

### Running Jupyter notebooks

To run a Jupyter notebook:

1. Start the Jupyter Notebook server:
    ```bash
    jupyter notebook
    ```
2. In your web browser, navigate to the provided local URL (usually `http://localhost:8888`).
3. Open the desired `.ipynb` notebook file and run the cells as needed.

____
