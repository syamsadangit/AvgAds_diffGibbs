import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#plt.style.use('_mpl-gallery')

def read_data_file(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', header=0)  
        print(df)

        return df
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

#READ FILES 
H1_ads = 'Diff-Ads-Energy-1H.dat'
H1_ads_df = read_data_file(H1_ads)
H1_ads_df.iloc[:, 0] = pd.to_numeric(H1_ads_df.iloc[:, 0], errors='coerce').astype(float)
H1_ads_df.iloc[:, 1] = pd.to_numeric(H1_ads_df.iloc[:, 1], errors='coerce').astype(float)

H2_ads = 'Diff-Ads-Energy-2H.dat'
H2_ads_df = read_data_file(H2_ads)
H2_ads_df.iloc[:, 0] = pd.to_numeric(H2_ads_df.iloc[:, 0], errors='coerce').astype(float)
H2_ads_df.iloc[:, 1] = pd.to_numeric(H2_ads_df.iloc[:, 1], errors='coerce').astype(float)

x  = H1_ads_df.iloc[:, 0]
y1 = H1_ads_df.iloc[:, 1]
y2 = H2_ads_df.iloc[:, 1]

print(H1_ads_df.dtypes)
print(H2_ads_df.dtypes)

print(y1, y2)

#PLOTTING PARAMETERS

bar_width = 0.2
frac = 0.5

# BAR POSITIONS FOR Y1 Y2
bar_positions_y1 = x #- (frac/3)
bar_positions_y2 = x #+ (frac/3)

# PLOTTING BAR PLOT
plt.bar(bar_positions_y2, y2, width=bar_width, color='orange', label='H2')
plt.bar(bar_positions_y1, y1, width=bar_width, color='blue', label='H1')

# Add labels and legend
plt.xlabel('Hydrogen Coverage')
plt.ylabel('Differential Adsorption Energy (eV)')
plt.legend()

# Show the plot
plt.savefig("H1-H2-Ads-comnpare.png")
plt.show()
