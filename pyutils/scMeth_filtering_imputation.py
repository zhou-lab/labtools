import pandas as pd
import numpy as np
import argparse
from collections import Counter
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA

def list_stats(my_list):
    # Get the counts of unique elements
    element_counts = Counter(my_list)

    # Print the unique elements and their counts
    for element, count in element_counts.items():
        print(f"Element: {element}, Count: {count}")

def filter_data(data, bin_num=10000, cell_pro=0.5):
    # Step 1: Filter cells based on non-NA counts
    data.columns = ['category'] + [f'feature_{i}' for i in range(1, data.shape[1])]
    non_na_counts = data.notna().sum(axis=1)
    filtered_data = data[non_na_counts >= bin_num]
    print("Finished cell filtering.")

    # Step 2: Filter columns based on threshold
    threshold = len(filtered_data) * cell_pro
    filtered_columns = filtered_data.columns[filtered_data.notna().sum() >= threshold]
    data_filtered = filtered_data[filtered_columns]
    data_filtered = data_filtered[['category'] + [col for col in filtered_columns if col != 'category']]
    print("Finished bin filtering.")

    # Step 3: Prepare data for joint imputation-clustering
    # Convert to numpy array for processing, excluding the 'category' column
    data_matrix = data_filtered.drop(columns=['category']).to_numpy()

    return data_filtered, data_matrix

def initial_imputation(data_filtered, strategy):
    # Step 4: Initial imputation (mean imputation)
    if strategy == "bin":
        imputer = SimpleImputer(strategy='mean')
        imputed_data = pd.DataFrame(imputer.fit_transform(data_filtered.iloc[:, 1:]))
    elif strategy == "mean":
        # Calculate the mean of all existing values (excluding the first column)
        overall_mean = data_filtered.iloc[:, 1:].stack().mean()
        # Impute missing values with this overall mean
        imputer = SimpleImputer(strategy='constant', fill_value=overall_mean)
        imputed_data = pd.DataFrame(imputer.fit_transform(data_filtered.iloc[:, 1:]))
    else:
        # Impute missing values (NAs) with a value
        imputer = SimpleImputer(strategy='constant', fill_value=float(strategy))
        imputed_data = pd.DataFrame(imputer.fit_transform(data_filtered.iloc[:, 1:]))
    return imputed_data
    
if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Joint imputation-clustering of DNA modification data.')
    parser.add_argument('input_file', type=str, help='Path to the input tab-separated file.')
    parser.add_argument('output_file', type=str, help='Path to the output tab-separated file.')
    parser.add_argument('bin_num', type=int, help='Number of bins a cell should have.')
    parser.add_argument('cell_proportion', type=float, help='Proportion of cells a bin should have.')
    parser.add_argument('initial_impute', type=str, help='The strategy for initial imputation.')
    args = parser.parse_args()

    data = pd.read_csv(args.input_file, sep='\t', header=None, low_memory=False)
    print("Finished loading tab-separated file.")
    data_filtered, data_matrix = filter_data(data, args.bin_num, args.cell_proportion)
    print("Finished filtering.")
    imputed_data = initial_imputation(data_filtered, args.initial_impute)

    imputed_df = pd.DataFrame(imputed_data)
    imputed_df['category'] = data_filtered['category'].values

    # Save to a CSV file
    imputed_df.to_csv(args.output_file, index=False)