import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# --- 1. THE SMART METADATA PARSER ---
def get_sample_groups(file_path):
    sample_mapping = {}
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        titles, geo_ids = [], []
        for line in lines:
            if line.startswith('!Sample_title') and 'aspartic' not in line.lower():
                titles = [t.strip('" \n') for t in line.split('\t')[1:]]
            if line.startswith('!Sample_geo_accession'):
                geo_ids = [g.strip('" \n') for g in line.split('\t')[1:]]

        if titles and geo_ids:
            for title, geo_id in zip(titles, geo_ids):
                clean_name = title.split(' ')[0].split('_')[0].split(':')[0].strip().capitalize()
                if clean_name not in sample_mapping:
                    sample_mapping[clean_name] = []
                sample_mapping[clean_name].append(geo_id)
        
        # Fallback if parsing fails to find Brain/Lung
        if "Brain" not in sample_mapping:
            sample_mapping = {
                "Brain": ["GSM15785", "GSM15786", "GSM15787", "GSM15788", "GSM15789"],
                "Lung": ["GSM15790", "GSM15791", "GSM15792", "GSM15794", "GSM15795"],
                "Colon": ["GSM15796", "GSM15797", "GSM15798", "GSM15799", "GSM15800"]
            }
        return sample_mapping
    except Exception as e:
        print(f"[!] Warning: Could not parse metadata. Using default categories. ({e})")
        return {"Brain": [], "Lung": []}

# --- 2. DATA LOADING ---
def load_and_clean_data(file_path):
    print(f"\n[1/2] Loading dataset: {file_path}...")
    try:
        # Check if file exists first
        df = pd.read_csv(file_path, sep='\t', comment='!', skiprows=30) 
        if df.empty:
            raise ValueError("The file appears to be empty.")
            
        df.rename(columns={df.columns[0]: 'Probe_ID'}, inplace=True)
        df['Probe_ID'] = df['Probe_ID'].astype(str).str.strip()
        df.set_index('Probe_ID', inplace=True)
        
        # Coerce non-numeric data to NaN and drop empty columns
        df = df.apply(pd.to_numeric, errors='coerce')
        df.dropna(how='all', axis=1, inplace=True)
        
        print(f"Success: Loaded {df.shape[0]} gene probes.")
        return df
    except FileNotFoundError:
        print(f"[!] Error: The file '{file_path}' was not found in this folder.")
    except pd.errors.EmptyDataError:
        print("[!] Error: No data found in the file.")
    except Exception as e:
        print(f"[!] An unexpected error occurred during loading: {e}")
    return None

# --- 3. ANALYSIS INTERFACE ---
def run_analysis(df, sample_map):
    try:
        print("\n" + "="*30 + "\n      ANALYSIS MENU\n" + "="*30)
        available_tissues = sorted(list(sample_map.keys()))
        print("Available Organs/Tissues:")
        print(" | ".join(available_tissues))
        print("-" * 30)
        
        target_probe = input("\nEnter Probe ID: ").strip()
        if target_probe not in df.index:
            print(f"[!] Error: Probe ID '{target_probe}' not found.")
            return

        g1 = input(f"Enter Group 1 Name: ").capitalize()
        g2 = input(f"Enter Group 2 Name: ").capitalize()

        if g1 not in sample_map or g2 not in sample_map:
            print("[!] Error: One of those organ names is not in the detected list.")
            return

        # Safely extract data
        data1 = df.loc[target_probe][sample_map[g1]].dropna()
        data2 = df.loc[target_probe][sample_map[g2]].dropna()

        if len(data1) < 2 or len(data2) < 2:
            print("[!] Error: Not enough data points in these groups to perform a T-test.")
            return

        t_stat, p_val = stats.ttest_ind(data1, data2, equal_var=False)
        
        print(f"\n--- Results for {target_probe} ---\n{g1} Mean: {data1.mean():.4f}\n{g2} Mean: {data2.mean():.4f}\nP-Value: {p_val:.6f}")

        # Visualization
        plt.figure(figsize=(8, 6))
        plot_df = pd.DataFrame({
            'Expression': pd.concat([data1, data2]),
            'Group': [g1]*len(data1) + [g2]*len(data2)
        })
        sns.boxplot(x='Group', y='Expression', data=plot_df, palette='Set2', hue='Group', legend=False)
        plt.title(f'Comparison: {target_probe}\n{g1}1 vs {g2} (p={p_val:.4f})')
        plt.show()

    except KeyError as e:
        print(f"[!] Error: Sample ID {e} was mentioned in metadata but is missing from data columns.")
    except Exception as e:
        print(f"[!] Error during analysis: {e}")

# --- 4. MAIN PROGRAM ---
def main():
    FILE_NAME = 'GSE1000_series_matrix.txt'
    data = load_and_clean_data(FILE_NAME)
    if data is None:
        input("\nPress Enter to exit...")
        return

    sample_map = get_sample_groups(FILE_NAME)

    while True:
        try:
            print("\n" + "*"*5 + "GENE EXPRESSION COMPARISION TOOL" + "*"*5)
            print("1. Search for a Gene/Probe ID\n2. Compare Two Organs\n3. Exit")
            choice = input("\nAction (1-3): ")
            
            if choice == '1':
                term = input("Search term: ").strip()
                matches = data.index[data.index.str.contains(term, case=False)].tolist()
                print(f"Top matches: ")
                for match in matches[:10]:
                    print(f" - {match}")
            elif choice == '2':
                run_analysis(data, sample_map)
            elif choice == '3':
                print("Thank you for using the Gene Expression Comparison Tool!")
                break
        except KeyboardInterrupt:
            print("\nProgram closed by user.")
            break
        except Exception as e:
            print(f"[!] Menu Error: {e}")

if __name__ == "__main__":
    main()
