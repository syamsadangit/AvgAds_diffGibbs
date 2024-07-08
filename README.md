# Overview

The Python code `Avg_diff_adsenergy.py` is designed to analyze vibrational energy, adsorption free energy, Gibbs free energy per hydrogen, and differential Gibbs free energy of hydrogen adsorption as a function of coverage on slab surfaces. The calculations and analysis are based on inputs of vibrational energy data, adsorption energy data, and other system-specific parameters.

# Usage

To use these functions, ensure you have the necessary input files containing vibrational energy and adsorption energy data.

## Necessary Inputs

### Vibrational Energy File

This file should contain the vibrational energy data needed for the `VibEneCollect` function. The file should be in a specific format where each line consists of vibrational energies followed by remarks, separated by spaces. The last element of each line should be the remark with a "#" character.

**Example format:**

0.123 0.456 0.789 #remark1



### Adsorption Energy File

This file should contain the adsorption energy data needed for the `AvgEads`, `diffGibbs`, `cleandf`, `write_AvgEads`, and `write_diffGibbs` functions. The file should be a CSV file with three columns: `nH`, `Eads`, and `Remarks`.

**Example format:**

nH,Eads,Remarks
1,0.123,remark1


### Reference Slab Energy

You will need the energy of the optimized adsorbed slab and the reference slab. These values are passed directly to the functions as variables, not as files. For example, `Eslab` and `EslabA7P` or `EslabAA7` are used in the `EneCorrectVib` function. These values should be available as part of your experimental or computational data.

## Ensure the following libraries are installed in your Python environment:

- numpy
- pandas
- matplotlib
- scipy

# Functions

### 1. `iseven(num)`

**Purpose:** Checks if a number is even.  
**Inputs:**  
- `num`: Integer to be checked.  
**Outputs:**  
- `condition`: Boolean indicating if `num` is even.  

### 2. `FreeEntVib(T, Evib)`

**Purpose:** Calculates the free energy of vibration, vibrational term, and entropy term based on temperature and vibrational energy.  
**Inputs:**  
- `T`: Temperature (K).  
- `Evib`: Vibrational energy (eV).  
**Outputs:**  
- `F`: Free energy of vibration.  
- `termvib`: Vibrational term.  
- `termentropy`: Entropy term.  

### 3. `write_vib_entropy(Vibenefile)`

**Purpose:** Writes the vibrational free energy and entropy due to vibration to a file.  
**Inputs:**  
- `Vibenefile`: Filename containing vibrational energy data.  
**Outputs:**  
- None (writes results to a file).  

### 4. `EneCorrectVib(Vibenefile, dfEads, Eslab)`

**Purpose:** Corrects the adsorption energy by adding the corresponding site-vibrational energy.  
**Inputs:**  
- `Vibenefile`: Filename containing vibrational energy data.  
- `dfEads`: DataFrame of adsorption energies.  
- `Eslab`: Energy of the reference slab.  
**Outputs:**  
- `dfeads`: DataFrame with corrected adsorption energies.  
- `Eslab`: Corrected slab energy.  

### 5. `VibEneCollect(input_file, dict)`

**Purpose:** Reads vibrational energy data from an input file into a DataFrame and dictionary.  
**Inputs:**  
- `input_file`: Filename containing vibrational energy data.  
- `dict`: Dictionary to store vibrational energy data.  
**Outputs:**  
- `df`: DataFrame containing vibrational energy and remarks.  
- `dict`: Updated dictionary with vibrational energy data.  

### 6. `AvgEads(df, Eslab, Vibenefile)`

**Purpose:** Calculates the average adsorption energy per hydrogen and Gibbs free energy per hydrogen.  
**Inputs:**  
- `df`: DataFrame of adsorption energies.  
- `Eslab`: Energy of the reference slab.  
- `Vibenefile`: Filename containing vibrational energy data.  
**Outputs:**  
- `df1`: DataFrame of average adsorption energies (no correction).  
- `df2`: DataFrame of Gibbs free energy per hydrogen (vibrational energy corrected).  

### 7. `diffGibbs(df1, df2, neibr)`

**Purpose:** Calculates the differential Gibbs free energy.  
**Inputs:**  
- `df1`: DataFrame of adsorption energies without reference energy values.  
- `df2`: DataFrame of adsorption energies with reference energy values.  
- `neibr`: Neighbor value (typically 1 or 2).  
**Outputs:**  
- `DGibbsdf`: DataFrame of differential Gibbs free energies.  

### 8. `cleandf(df, Eslab, Vibenefile, rem_avail=False)`

**Purpose:** Cleans the DataFrame to find the minimum energy adsorption configuration for a given number of hydrogen atoms.  
**Inputs:**  
- `df`: DataFrame of adsorption energies.  
- `Eslab`: Energy of the reference slab.  
- `Vibenefile`: Filename containing vibrational energy data.  
- `rem_avail`: Boolean to include/removal of specific remarks.  
**Outputs:**  
- `dfdr1`: Cleaned DataFrame with unique configurations.  
- `dfclean`: Cleaned DataFrame including reference energy.  

### 9. `write_AvgEads(input_file, Eslab, system, Vibenefile)`

**Purpose:** Writes the average adsorption energy and Gibbs free energy per hydrogen to files.  
**Inputs:**  
- `input_file`: Filename containing adsorption energy data.  
- `Eslab`: Energy of the reference slab.  
- `system`: System name/identifier.  
- `Vibenefile`: Filename containing vibrational energy data.  
**Outputs:**  
- None (writes results to files).  

### 10. `write_diffGibbs(input_file, Eslab, system, Vibenefile)`

**Purpose:** Writes the differential Gibbs free energy to files.  
**Inputs:**  
- `input_file`: Filename containing adsorption energy data.  
- `Eslab`: Energy of the reference slab.  
- `system`: System name/identifier.  
- `Vibenefile`: Filename containing vibrational energy data.  
**Outputs:**  
- None (writes results to files).  

### 11. `plotAvgAdsene(dfEads, df, input_file, system)`

**Purpose:** Plots the average adsorption energy and Gibbs free energy per hydrogen.  
**Inputs:**  
- `dfEads`: DataFrame of average adsorption energies.  
- `df`: DataFrame of Gibbs free energy per hydrogen.  
- `input_file`: Filename containing adsorption energy data.  
- `system`: System name/identifier.  
**Outputs:**  
- None (saves plot to a file).  

### 12. `plotdGibbs(df, input_file, system)`

**Purpose:** Plots the differential Gibbs free energy.  
**Inputs:**  
- `df`: DataFrame of differential Gibbs free energies.  
- `input_file`: Filename containing adsorption energy data.  
- `system`: System name/identifier.  
**Outputs:**  
- None (saves plot to a file).  

### 13. `plotdAdsEne(df, input_file, system)`

**Purpose:** Plots the differential adsorption energy.  
**Inputs:**  
- `df`: DataFrame of differential adsorption energies.  
- `input_file`: Filename containing adsorption energy data.  
- `system`: System name/identifier.  
**Outputs:**  
- None (saves plot to a file).

