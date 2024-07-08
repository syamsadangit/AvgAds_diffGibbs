import pandas as pd

def remove_hash_from_last_column(file_name):
    df = pd.read_csv(file_name, delimiter='\s+')  # Assuming the files are space-delimited
    df.iloc[:, -1] = df.iloc[:, -1].apply(lambda x: x.replace('#', ''))
    return df

file_names = ["AA7cvgGibbsFene.dat", "AA7+PcvgGibbsFene.dat", "GibbsFeneminAA7cvg.dat", "GibbsFeneminAA7+Pcvg.dat"]

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
    num_columns = len(df.columns) - 1
    headers = 'No. of H (n) & $\Delta$G(n) & Sites of adsorption'

    rows = ""
    for _, row in df.iterrows():
        row_data = [str(row.get(col, '')) for idx, col in enumerate(df.columns) if idx != 1]
        rows += ' & '.join(row_data) + r' \\ \hline'
        rows += '\n'

    if len(df) > 30:  # Number of rows to decide whether to use longtable or regular table
        return template_longtable.replace('<column_specification>', '|'.join(['l'] * num_columns))\
                                 .replace('<caption>', f'A sample long table (Table {table_number})')\
                                 .replace('<label>', f'tab:long{table_number}')\
                                 .replace('<num_columns>', str(num_columns))\
                                 .replace('<headers>', headers)\
                                 .replace('<rows>', rows)
    else:
        return template_table.replace('<column_specification>', '|'.join(['l'] * num_columns))\
                             .replace('<caption>', f'A sample table (Table {table_number})')\
                             .replace('<label>', f'tab:table{table_number}')\
                             .replace('<headers>', headers)\
                             .replace('<rows>', rows)

for idx, file_name in enumerate(file_names, start=1):
    df = remove_hash_from_last_column(file_name)
    output_file_name = file_name.replace(".dat", ".tex")
    latex_table = generate_longtable_or_table(df, table_number=idx)

    # Opening the output .tex file
    with open(output_file_name, 'w') as f:
        f.write(latex_table)

