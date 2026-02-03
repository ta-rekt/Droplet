"""
Visualization of Flory-Huggins Prewetting Transition
Heatmap of largest_droplet_size averaged over final 20% of frames
for varying temperature T and concentration α
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import sys

# Configuration
DATA_DIR = Path("/home/tohme/.julia/dev/Droplet/data")
TEMPERATURES = [0.56, 0.58, 0.60, 0.62, 0.64, 0.66]
NCPU = 32  # Number of CPU runs per temperature (32 concentrations)

def load_and_aggregate_data(date_folders):
    """
    Load JSON data from all date/temperature combinations and aggregate.
    
    Parameters:
        date_folders: list of date folder names to aggregate
    
    Returns:
        heatmap_data: (n_temps, n_cpus) array of averaged largest_droplet_size
        concentrations: list of concentration values from filenames
    """
    
    # Dictionary to store aggregated data: {temp: {cpu_idx: [values]}}
    aggregated = {t: {i: [] for i in range(1, NCPU + 1)} for t in TEMPERATURES}
    
    print(f"Scanning data directories: {date_folders}")
    print(f"Temperatures: {TEMPERATURES}")
    print(f"Expected CPU runs per temperature: {NCPU}")
    
    count_by_temp = {t: 0 for t in TEMPERATURES}
    
    # Iterate through date folders
    for date_folder in date_folders:
        date_path = DATA_DIR / date_folder
        if not date_path.exists():
            print(f"Warning: {date_path} does not exist")
            continue
            
        # Iterate through temperature folders
        for temp in TEMPERATURES:
            temp_folder_pattern = f"flory_huggins_prewetting_42_T_{temp:.2f}"
            temp_path = date_path / temp_folder_pattern
            
            if not temp_path.exists():
                print(f"Warning: {temp_path} does not exist")
                continue
            
            # Load data from each CPU output file
            for cpu_idx in range(1, NCPU + 1):
                json_file = temp_path / f"output_cpu{cpu_idx}.json"
                
                if not json_file.exists():
                    print(f"Warning: {json_file} does not exist")
                    continue
                
                try:
                    with open(json_file, 'r') as f:
                        data = json.load(f)
                    
                    # Data is nested with two numeric keys at top level
                    first_key = list(data.keys())[0]
                    second_key = list(data[first_key].keys())[0]
                    sim_data = data[first_key][second_key]
                    
                    # Extract necessary fields
                    largest_droplet_size = sim_data['largest_droplet_size']
                    nframes = sim_data['nframes']
                    concentration = sim_data['concentration']
                    
                    # Calculate threshold for final 20% of frames
                    threshold_frame = int(0.8 * nframes)
                    
                    # Average over final 20% of frames
                    avg_size = np.mean(largest_droplet_size[threshold_frame:])
                    
                    # Store aggregated value
                    aggregated[temp][cpu_idx].append({
                        'avg_size': avg_size,
                        'concentration': concentration,
                        'date': date_folder
                    })
                    
                    count_by_temp[temp] += 1
                    
                except Exception as e:
                    print(f"Error processing {json_file}: {e}")
                    continue
    
    print(f"\nData loading summary:")
    for temp, count in count_by_temp.items():
        print(f"  T={temp}: {count} files loaded")
    
    # Now organize data into heatmap format
    # Shape: (n_temps, n_cpus) where each element is the average over dates
    heatmap_data = np.zeros((len(TEMPERATURES), NCPU))
    concentrations_by_cpu = {}
    
    for temp_idx, temp in enumerate(TEMPERATURES):
        for cpu_idx in range(1, NCPU + 1):
            values = aggregated[temp][cpu_idx]
            if values:
                # Average across dates
                avg_across_dates = np.mean([v['avg_size'] for v in values])
                heatmap_data[temp_idx, cpu_idx - 1] = np.log(avg_across_dates)
                
                # Store concentration (should be same for all dates)
                if cpu_idx not in concentrations_by_cpu:
                    concentrations_by_cpu[cpu_idx] = values[0]['concentration']
            else:
                heatmap_data[temp_idx, cpu_idx - 1] = np.nan
    
    # Extract concentrations in order
    concentrations = [concentrations_by_cpu[i] for i in range(1, NCPU + 1)]
    
    # Sort by concentration to ensure proper ordering
    sorted_indices = np.argsort(concentrations)
    heatmap_data_sorted = heatmap_data[:, sorted_indices]
    concentrations_sorted = np.array(concentrations)[sorted_indices]
    
    return heatmap_data_sorted, concentrations_sorted

def plot_heatmap(heatmap_data, concentrations, temperatures, date_str, file_name):
    """
    Create and save heatmap visualization.
    
    Parameters:
        heatmap_data: (n_temps, n_concs) array
        concentrations: array of concentration values
        temperatures: list of temperature values
        date_str: date folder name to append to filename
    """
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Create heatmap
    im = ax.imshow(heatmap_data, aspect='auto', cmap='viridis', origin='lower',
                   extent=[concentrations[0], concentrations[-1], 
                          temperatures[0], temperatures[-1]],
                   interpolation='nearest', vmin=0, vmax=7)
    
    # Labels and title
    ax.set_xlabel('Concentration (α)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Temperature (T)', fontsize=12, fontweight='bold')
    ax.set_title('Flory-Huggins Prewetting Transition\nAverage Largest Droplet Size (Final 20% of Frames)',
                fontsize=13, fontweight='bold', pad=20)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label='Largest Droplet Size (lattice units)')
    
    # Set tick positions
    ax.set_yticks(temperatures)
    ax.set_xticks(np.linspace(concentrations[0], concentrations[-1], 5))
    
    # Formatting
    ax.grid(False)
    fig.tight_layout()
    
    # Save figure
    output_path = Path(f"/home/tohme/.julia/dev/Droplet/media/figures/{date_str}/{file_name}.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nHeatmap saved to: {output_path}")
    
    plt.close()

def main():
    print("=" * 70)
    print("Flory-Huggins Prewetting Transition Visualization")
    print("=" * 70)
    
    # Get date folder from command-line argument
    if len(sys.argv) < 2:
        print("Usage: python prewetting_heatmap.py <date_folder>")
        print("Example: python prewetting_heatmap.py 20260131")
        sys.exit(1)
    
    date_folders = [sys.argv[1]]
    file_name = sys.argv[2]
    
    # Load and aggregate data
    heatmap_data, concentrations = load_and_aggregate_data(date_folders)
    
    # Print statistics
    print(f"\nData statistics:")
    print(f"  Heatmap shape: {heatmap_data.shape}")
    print(f"  Temperature range: {TEMPERATURES[0]:.2f} - {TEMPERATURES[-1]:.2f}")
    print(f"  Concentration range: {concentrations[0]:.6f} - {concentrations[-1]:.6f}")
    print(f"  Mean droplet size: {np.nanmean(heatmap_data):.2f}")
    print(f"  Std droplet size: {np.nanstd(heatmap_data):.2f}")
    print(f"  Min droplet size: {np.nanmin(heatmap_data):.2f}")
    print(f"  Max droplet size: {np.nanmax(heatmap_data):.2f}")
    
    # Create and save heatmap
    plot_heatmap(heatmap_data, concentrations, TEMPERATURES, date_folders[0], file_name)
    
    print("\n" + "=" * 70)
    print("Visualization complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
