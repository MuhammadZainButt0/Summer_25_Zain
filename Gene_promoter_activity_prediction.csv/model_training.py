import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import matplotlib.pyplot as plt

def train_evaluate_model(input_file="data/Final_dataset.csv"):
    """
    Trains an XGBoost model to predict gene expression (FPKM) and creates visualizations.

    Args:
        input_file (str): The path to the final dataset for training.

    Returns:
        pandas.DataFrame or None: A DataFrame with model predictions, or None if an error occurs.
    """
    # Check if the input file exists. If not, print an error and exit.
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file '{input_file}' not found. Please ensure the previous pipeline steps completed successfully.")
        return None

    # Load the dataset and remove any rows with missing values.
    df = pd.read_csv(input_file)
    df = df.dropna()

    # Identify columns that are not features for the model.
    non_feature_cols = ['gene_id', 'promoter_sequence', 'expected_count', 'TPM', 'FPKM']
    # Create a list of feature columns.
    feature_cols = [col for col in df.columns if col not in non_feature_cols]

    # Identify and exclude any non-numeric columns from the features.
    non_numeric_cols = [col for col in feature_cols if not pd.api.types.is_numeric_dtype(df[col])]
    if non_numeric_cols:
        print(f"⚠️ Warning: Non-numeric columns found: {non_numeric_cols}. These will be excluded.")
        feature_cols = [col for col in feature_cols if col not in non_numeric_cols]

    # If no feature columns are left, print an error and exit.
    if not feature_cols:
        print("❌ Error: No numeric feature columns available for training.")
        return None

    # Define the features (X) and the target variable (y).
    X = df[feature_cols]
    y = np.log1p(df['FPKM'])  # Use a log transformation on FPKM for better model performance.

    # Select the top 50 features that are most correlated with the target variable.
    try:
        selector = SelectKBest(score_func=f_regression, k=min(50, len(feature_cols)))
        X_selected = selector.fit_transform(X, y)
        selected_features = X.columns[selector.get_support()].tolist()
    except Exception as e:
        print(f"❌ Error during feature selection: {e}")
        return None

    # Scale the feature values to have a mean of 0 and a standard deviation of 1.
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_selected)

    # Split the data into training and testing sets.
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

    # Set up the XGBoost model and a grid of parameters for tuning.
    xgb = XGBRegressor(random_state=42)
    param_grid = {
        'n_estimators': [100, 200],
        'max_depth': [3, 5, 7],
        'learning_rate': [0.01, 0.1],
        'subsample': [0.8, 1.0]
    }
    # Use GridSearchCV to find the best model parameters based on R-squared scoring.
    grid_search = GridSearchCV(xgb, param_grid, cv=5, scoring='r2', n_jobs=-1)
    grid_search.fit(X_train, y_train)

    # Get the best model from the grid search and print its parameters.
    best_model = grid_search.best_estimator_
    print("Best Parameters:", grid_search.best_params_)
    
    # Make predictions on the test set and evaluate the model's performance.
    y_pred = best_model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    print(f"Mean Squared Error: {mse:.4f}")
    print(f"Root Mean Squared Error: {rmse:.4f}")
    print(f"Mean Absolute Error: {mae:.4f}")
    print(f"R² Score: {r2:.4f}")

    # Plot the top 10 most important features for the model's predictions.
    importances = best_model.feature_importances_
    feature_importance = pd.Series(importances, index=selected_features).sort_values(ascending=False)
    plt.figure(figsize=(10, 6))
    feature_importance.head(10).plot(kind='bar')
    plt.title('Top 10 Feature Importances for Predicting FPKM')
    plt.xlabel('Features')
    plt.ylabel('Importance')
    plt.tight_layout()
    plt.savefig('data/feature_importances.png')
    plt.close()

    # Create a scatter plot to compare the actual and predicted FPKM values.
    results = pd.DataFrame({'Actual_FPKM': np.expm1(y_test), 'Predicted_FPKM': np.expm1(y_pred)})
    plt.figure(figsize=(8, 6))
    plt.scatter(results['Actual_FPKM'], results['Predicted_FPKM'], alpha=0.5)
    plt.plot([0, results['Actual_FPKM'].max()], [0, results['Actual_FPKM'].max()], 'r--')
    plt.xlabel('Actual FPKM')
    plt.ylabel('Predicted FPKM')
    plt.title('Predicted vs Actual FPKM')
    plt.grid(True)
    plt.savefig('data/fpkm_pred_vs_actual.png')
    plt.close()

    # Save the prediction results to a CSV file.
    results.to_csv('data/predictions.csv', index=False)
    print("✅ Predictions saved to 'data/predictions.csv'")
    
    return results

def delete_csvs_except_final(data_dir="data"):
    """
    Deletes all CSV files in the data directory except for a few specified files.

    Args:
        data_dir (str): The directory where the data files are stored.
    """
    # Check if the data directory exists.
    if not os.path.exists(data_dir):
        print(f"❌ Error: Directory '{data_dir}' not found.")
        return

    # A list of files to keep.
    keep_files = ["Final_dataset.csv", "predictions.csv", "promoter_sequences.csv"]

    # Loop through all files in the directory.
    for file in os.listdir(data_dir):
        # If a file is a CSV and not in the 'keep_files' list, delete it.
        if file.endswith(".csv") and file not in keep_files:
            file_path = os.path.join(data_dir, file)
            try:
                os.remove(file_path)
            except Exception as e:
                # If there's an error during deletion, do nothing and continue.
                pass

# This allows the delete_csvs_except_final function to be run directly from the command line.
if __name__ == "__main__":
    delete_csvs_except_final() 