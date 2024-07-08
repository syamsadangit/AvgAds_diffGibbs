import pandas as pd
import numpy as np

def remove_hash_from_last_column(file_name):
    """
    Reads a space-delimited file into a DataFrame and removes any '#' characters 
    from the last column.
    
    Parameters:
        file_name (str): The name of the file to read.
        
    Returns:
        pd.DataFrame: The cleaned DataFrame.
    """
    df = pd.read_csv(file_name, delimiter='\s+')
    df.iloc[:, -1] = df.iloc[:, -1].apply(lambda x: x.replace('#', ''))
    return df


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
    df.iloc[:, 2] = pd.Series(np.round(df.iloc[:, 2].to_numpy(), 2))

    # Sort the DataFrame first by the first column and then by the second column
    df_sorted = df.sort_values(by=[df.columns[0], df.columns[2]])
    

    return df_sorted


# List of filenames to process
file_names = ["AA7cvgGibbsFene.dat", "AA7+PcvgGibbsFene.dat", "GibbsFeneminAA7cvg.dat", "GibbsFeneminAA7+Pcvg.dat"]

# LaTeX table templates
template_longtable = r'''
\clearpage
\begin{center}
\centering
\begin{longtable}{|<column_specification>|}
\caption{<caption>}
\label{<label>} \\
\hline
<headers> \\ \hline
\endfirsthead

\multicolumn{<num_columns>}{c}{\tablename\ \thetable{} -- Continued from previous page} \\ \hline
<headers> \\ \hline
\endhead

\hline
\endfoot

\hline
\endlastfoot
<rows>
\end{longtable}
\end{center}
'''

template_table = r'''
\begin{table}[ht]
\centering
\caption{<caption>}
\label{<label>}
\begin{tabular}{|<column_specification>|}
\hline
<headers> \\ \hline
<rows>
\end{tabular}
\end{table}
'''

def generate_longtable_or_table(df, table_number):
    """
    Generates a LaTeX table or longtable from a DataFrame.
    
    Parameters:
        df (pd.DataFrame): The DataFrame to convert to a LaTeX table.
        table_number (int): The table number to use in the LaTeX caption and label.
        
    Returns:
        str: The LaTeX table as a string.
    """
    num_columns = len(df.columns) - 1
    headers = ' & '.join(df.columns.drop(df.columns[1]))

    rows = ""
    for _, row in df.iterrows():
        row_data = [str(row.get(col, '')) for idx, col in enumerate(df.columns) if idx != 1]
        rows += ' & '.join(row_data) + r' \\ \hline'
        rows += '\n'

    if len(df) > 30:  # Use longtable if there are more than 30 rows
        return template_longtable.replace('<column_specification>', '|'.join(['l'] * num_columns))\
                                 .replace('<caption>', f'A sample long table (Table {table_number})')\
                                 .replace('<label>', f'tab:long{table_number}')\
                                 .replace('<num_columns>', str(num_columns))\
                                 .replace('<headers>', headers)\
                                 .replace('<rows>', rows)
    else:  # Use regular table for fewer than 30 rows
        return template_table.replace('<column_specification>', '|'.join(['l'] * num_columns))\
                             .replace('<caption>', f'A sample table (Table {table_number})')\
                             .replace('<label>', f'tab:table{table_number}')\
                             .replace('<headers>', headers)\
                             .replace('<rows>', rows)

# Process each file and generate the corresponding LaTeX table
for idx, file_name in enumerate(file_names, start=1):
    df = remove_hash_from_last_column(file_name)
    df_sorted_rounded = sort_and_round_dataframe(df)
    output_file_name = file_name.replace(".dat", ".tex")
    latex_table = generate_longtable_or_table(df_sorted_rounded, table_number=idx)

    # Write the LaTeX table to a .tex file
    with open(output_file_name, 'w') as f:
        f.write(latex_table)

# Print brief documentation after running the script
print("""
Purpose of the Code:
This script processes space-delimited data files, removes any '#' characters from the last column, sorts the data based on the first column and then by the second column, rounds the second column to two decimal places, and generates LaTeX tables based on the cleaned and sorted data. The output LaTeX tables are written to new .tex files.

Output Data:
The following .tex files are generated:
- AA7cvgGibbsFene.tex
- AA7+PcvgGibbsFene.tex
- GibbsFeneminAA7cvg.tex
- GibbsFeneminAA7+Pcvg.tex
""")

