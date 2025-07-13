# Gene Expression Prediction Pipeline

This project contains simplified Python scripts for predicting gene expression levels based on promoter sequences. The code is designed to be beginner-friendly with clear comments explaining each step.

## Files Overview

- **`main.py`**: The main script that runs the entire pipeline
- **`data_preparation.py`**: Handles data conversion and merging
- **`sequence_fetcher.py`**: Fetches promoter sequences from Ensembl database
- **`feature_extraction.py`**: Extracts various features from DNA sequences
- **`one_hot_encoding.py`**: Converts DNA sequences to numerical format
- **`model_training.py`**: Trains and evaluates the machine learning model

## Setup Instructions

1. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Prepare your data:**
   - Place your gene expression TSV file in the `data/` folder
   - Name it `gene_expression.A549.biorep1.tsv`
   - The file should have columns: `gene_id`, `expected_count`, `TPM`, `FPKM`

3. **Run the pipeline:**
   ```bash
   python main.py
   ```

## What Each Script Does

### 1. Data Preparation (`data_preparation.py`)
- Converts TSV files to CSV format
- Cleans up gene IDs by removing version numbers
- Merges expression data with promoter sequences
- Creates the final dataset for training

### 2. Sequence Fetcher (`sequence_fetcher.py`)
- Fetches promoter sequences from the Ensembl database
- Handles API requests with proper error handling
- Saves sequences to CSV format

### 3. Feature Extraction (`feature_extraction.py`)
- Calculates GC content, AT content, and CpG density
- Computes Shannon entropy for sequence complexity
- Extracts k-mer frequencies (3-mers)
- Searches for regulatory motifs (TATA box, CAAT box, etc.)

### 4. One-Hot Encoding (`one_hot_encoding.py`)
- Converts DNA sequences to numerical format
- Each base (A, C, G, T) gets a unique binary code
- Creates fixed-length representations for machine learning

### 5. Model Training (`model_training.py`)
- Uses XGBoost for regression
- Performs feature selection
- Creates visualizations of results
- Saves predictions to CSV

## Output Files

After running the pipeline, you'll find these files in the `data/` folder:
- `Final_dataset.csv`: The complete dataset for training
- `predictions.csv`: Model predictions
- `feature_importances.png`: Visualization of important features
- `fpkm_pred_vs_actual.png`: Comparison of predicted vs actual values

## Key Features for Beginners

- **Clear Comments**: Every function and major step is explained
- **Error Handling**: Scripts check for missing files and handle errors gracefully
- **Progress Indicators**: Uses emojis and clear messages to show progress
- **Modular Design**: Each script can be run independently
- **Default Parameters**: Sensible defaults for all functions

## Troubleshooting

- **Missing input file**: Make sure your TSV file is in the `data/` folder
- **Import errors**: Install all dependencies with `pip install -r requirements.txt`
- **API errors**: The sequence fetcher includes delays to avoid overwhelming the server
- **Memory issues**: The pipeline processes data in chunks to manage memory usage

## Customization

You can modify the scripts by:
- Changing the list of gene IDs in `main.py`
- Adjusting feature extraction parameters in `feature_extraction.py`
- Modifying model parameters in `model_training.py`
- Adding new features or motifs as needed 