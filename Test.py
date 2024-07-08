import numpy as np
import pandas as pd

def sort_and_round_dataframe(df):
    """
    Sorts the DataFrame first by the first column and then by the second column,
    and rounds the second column to two decimal places.

    Parameters:
        df (pd.DataFrame): The DataFrame to sort and round.

    Returns:
        pd.DataFrame: The sorted and rounded DataFrame.
    """
    # Access the second column and round its values to two decimal places
    df.iloc[:, 1] = pd.Series(np.round(df.iloc[:, 1].to_numpy(), 2))
    
    # Sort the DataFrame first by the first column and then by the second column
    df_sorted = df.sort_values(by=[df.columns[0], df.columns[1]])
    
    return df_sorted

# Example usage:
data = {
    'Column1': [3, 1, 2],
    'Column2': [-0.12345, -0.67890, -0.23456],
    'Column3': [7, 8, 9]
}
df = pd.DataFrame(data)

# Before sorting and rounding
print("Before sorting and rounding:")
print(df)

# Apply sorting and rounding
df_sorted_rounded = sort_and_round_dataframe(df)

# After sorting and rounding
print("After sorting and rounding:")
print(df_sorted_rounded)

