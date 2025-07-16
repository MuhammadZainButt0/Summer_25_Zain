
####   *Methodlogy*



## 1. Project Overview

This project implements a deep learning pipeline to classify single nucleotide variants (SNVs) in human genes as *deleterious* or *neutral* based on their functional impact. The pipeline is modularized for clarity, maintainability, and reusability. The goal is to predict the class label of each SNP variant based on both sequence and engineered features.

---

## 2. Project Structure

The codebase is organized into the following modules:

- **`data_preprocessing.py`**  
  Handles data loading, feature engineering, and preprocessing.  
  - Reads the dataset (`Update.csv`)
  - Encodes categorical features (e.g., Chromosome, Ref, Alt, Variant_Type)
  - Computes and includes sequence-based features (Shannon Entropy, GC/AT content)
  - One-hot encodes the DNA sequence context
  - Returns a clean, numeric feature matrix and label vector

- **`dataset.py`**  
  Contains the `DNADataset` class, a PyTorch Dataset for efficient data loading.

- **`model.py`**  
  Defines the `DNAClassifier` neural network, a simple feedforward model for classification.

- **`train_eval.py`**  
  Contains functions for training (`train_model`), evaluation (`evaluate_model`), and visualization (`plot_confusion_matrix`).

- **`main.py`**  
  The entry point. Orchestrates the full pipeline: data loading, splitting, model training, evaluation, and reporting.

---
##3. Data Acquisition and Preprocessing

1. Data Source
To ensure the use of real, up-to-date, and publicly available data, we retrieved single nucleotide variant (SNV) information from the ClinVar database. Specifically, we downloaded the latest human SNP VCF file for the GRCh38 genome build from the official ClinVar FTP site.
ClinVar VCF Download:
The VCF file contains annotated SNVs, including their genomic positions, reference and alternate alleles, and clinical significance.
Reference Genome:
The human reference genome (GRCh38) in FASTA format was used to extract the sequence context around each SNV. (Due to its large size, this file was downloaded manually from the NCBI or Ensembl FTP site.)

2. SNP Mapping and Context Extraction
We parsed the ClinVar VCF file to identify  SNPs (variants where both the reference and alternate alleles are single nucleotides).
For each SNP, we mapped its position to the GRCh38 reference genome.
We extracted a 10 base pair (bp) context sequence for each SNP, consisting of 5 bp upstream, the SNP itself, and 4 bp downstream. This context provides local sequence information that may influence the functional impact of the variant.

3. Feature Engineering
For each extracted context sequence, we computed the following features:
Shannon Entropy: Quantifies the sequence complexity and diversity of nucleotides.
GC Content (%): The percentage of guanine (G) and cytosine (C) bases in the context.
AT Content (%): The percentage of adenine (A) and thymine (T) bases in the context.
These features, along with the variant annotations (chromosome, position, reference/alternate alleles), were compiled into a structured dataset for downstream machine learning.

4. Data Output
The final processed dataset was saved as a CSV file(Update.csv), which serves as the input for the deep learning classification pipeline.


## 4. Pipeline Steps

1. **Data Loading & Preprocessing**
   - The dataset is loaded from `Update.csv`.
   - Features are engineered and encoded as described above.
   - The final feature matrix is fully numeric and suitable for machine learning.

2. **Dataset & DataLoader**
   - The data is split into training and testing sets.
   - PyTorch `Dataset` and `DataLoader` objects are created for efficient batching.

3. **Model Definition**
   - A feedforward neural network is defined with one hidden layer and dropout for regularization.

4. **Training**
   - The model is trained using cross-entropy loss and the Adam optimizer.
   - Training progress (loss and accuracy) is printed every few epochs.

5. **Evaluation**
   - The trained model is evaluated on the test set.
   - Metrics reported: accuracy, precision, recall, F1-score, and a confusion matrix.

6. **Visualization**
   - The confusion matrix is plotted for visual inspection of classification performance.

---

## 5. How to Run

1. **Install dependencies** (if not already installed):
   ```bash
   pip install torch pandas numpy scikit-learn matplotlib seaborn
   ```

2. **Run the pipeline:**
   ```bash
   python main.py
   ```

---

## 6. *Findings*

- **Modularization:**  
  Breaking the code into modules improves readability, maintainability, and reusability.

- **Feature Engineering:**  
  Combining biological sequence features with categorical and numeric data can improve model performance.

- **PyTorch Workflow:**  
  Using PyTorch’s Dataset and DataLoader classes streamlines the training process for large datasets.

- **Evaluation:**  
  Multiple metrics(Classification report and accuracy) and visualizations (like the confusion matrix) provide a comprehensive view of model performance.


---

**Prepared by:**  
[Muhammad Zain Butt]  
[16-07-2025] 