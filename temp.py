import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.ndimage import uniform_filter1d


def compute_effective_force(dirname, keya_values, keyb=1, 
                            burn_in_frac=0.5, n_samples=10, thinning=100,
                            smoothing_window=3):
    """
    Compute effective force between prewet filaments from Hamiltonian vs distance.
    
    Parameters
    ----------
    dirname : str
        Directory containing simulation output files
    keya_values : list or array
        List of keya indices (corresponding to different distances)
    keyb : int
        Beta index (usually 1)
    burn_in_frac : float
        Fraction of frames to discard as burn-in
    n_samples : int
        Number of equilibrated samples to take per simulation
    thinning : int
        Number of frames between samples
    smoothing_window : int
        Window size for smoothing (must be odd)
        
    Returns
    -------
    distances : array
        Inter-filament distances
    mean_hams : array
        Mean Hamiltonian at each distance
    std_hams : array
        Standard deviation of Hamiltonian
    ham_smooth : array
        Smoothed Hamiltonian values
    forces : array
        Effective force (negative derivative of energy)
    """
    distances = []
    mean_hams = []
    std_hams = []
    
    for keya in keya_values:
        try:
            res_or_path = 'data/{}/output_cpu{}.json'.format(dirname, keya)
            with open(res_or_path, 'r') as f:
                res_ = json.load(f)
                res = res_[str(keya)][str(keyb)]
            
            ham = np.array(res.get("hamiltonian", []))
            distance = res.get("beta", keya)  # beta stores the distance
            
            if len(ham) == 0:
                continue
                
            # Compute burn-in index
            burn_in_idx = int(len(ham) * burn_in_frac)
            
            # Sample equilibrated values with thinning
            equilibrated = ham[burn_in_idx:]
            
            # Take n_samples evenly spaced (with thinning)
            if len(equilibrated) > n_samples * thinning:
                sample_indices = np.arange(0, n_samples * thinning, thinning)
                samples = equilibrated[sample_indices]
            else:
                # If not enough data, take what we can
                samples = equilibrated[::max(1, len(equilibrated) // n_samples)]
            
            distances.append(distance)
            mean_hams.append(np.mean(samples))
            std_hams.append(np.std(samples) / np.sqrt(len(samples)))  # standard error
            
        except Exception as e:
            print(f"Error loading keya={keya}: {e}")
            continue
    
    # Convert to arrays and sort by distance
    distances = np.array(distances)
    mean_hams = np.array(mean_hams)
    std_hams = np.array(std_hams)
    
    sort_idx = np.argsort(distances)
    distances = distances[sort_idx]
    mean_hams = mean_hams[sort_idx]
    std_hams = std_hams[sort_idx]
    
    # Smooth the Hamiltonian curve
    if smoothing_window > 1 and len(mean_hams) >= smoothing_window:
        ham_smooth = uniform_filter1d(mean_hams, size=smoothing_window, mode='nearest')
    else:
        ham_smooth = mean_hams.copy()
    
    # Compute force as negative derivative: F = -dH/dr
    # Use central differences where possible
    forces = np.zeros_like(ham_smooth)
    if len(distances) > 2:
        # Central difference for interior points
        for i in range(1, len(distances) - 1):
            dr = distances[i+1] - distances[i-1]
            if dr > 0:
                forces[i] = -(ham_smooth[i+1] - ham_smooth[i-1]) / dr
        
        # Forward/backward difference for endpoints
        if distances[1] != distances[0]:
            forces[0] = -(ham_smooth[1] - ham_smooth[0]) / (distances[1] - distances[0])
        if distances[-1] != distances[-2]:
            forces[-1] = -(ham_smooth[-1] - ham_smooth[-2]) / (distances[-1] - distances[-2])
    
    return distances, mean_hams, std_hams, ham_smooth, forces


def plot_effective_force(dirname, keya_values, keyb=1, 
                         burn_in_frac=0.5, n_samples=10, thinning=100,
                         smoothing_window=3, figsize=(12, 5)):
    """
    Plot Hamiltonian vs distance and effective force between filaments.
    """
    distances, mean_hams, std_hams, ham_smooth, forces = compute_effective_force(
        dirname, keya_values, keyb, burn_in_frac, n_samples, thinning, smoothing_window
    )
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Plot 1: Raw Hamiltonian vs distance
    axes[0].errorbar(distances, mean_hams, yerr=std_hams, fmt='o', 
                     color='C0', capsize=3, markersize=5, label='Mean H')
    axes[0].plot(distances, ham_smooth, '-', color='C1', linewidth=2, 
                 label=f'Smoothed (w={smoothing_window})')
    axes[0].set_xlabel('Inter-filament distance')
    axes[0].set_ylabel('Hamiltonian')
    axes[0].set_title('Energy vs Distance')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Effective force vs distance
    axes[1].plot(distances, forces, 'o-', color='C2', markersize=5, linewidth=2)
    axes[1].axhline(y=0, color='black', linestyle='--', alpha=0.5)
    axes[1].set_xlabel('Inter-filament distance')
    axes[1].set_ylabel('Force (F = -dH/dr)')
    axes[1].set_title('Effective Inter-filament Force')
    axes[1].grid(True, alpha=0.3)
    
    # Annotate attractive/repulsive regions
    if np.any(forces > 0) and np.any(forces < 0):
        axes[1].fill_between(distances, forces, 0, where=(forces > 0), 
                            alpha=0.3, color='red', label='Repulsive')
        axes[1].fill_between(distances, forces, 0, where=(forces < 0), 
                            alpha=0.3, color='blue', label='Attractive')
        axes[1].legend()
    
    # Plot 3: Potential well visualization
    # Shift energy so minimum is at zero
    ham_shifted = ham_smooth - np.min(ham_smooth)
    axes[2].plot(distances, ham_shifted, 'o-', color='C3', markersize=5, linewidth=2)
    axes[2].set_xlabel('Inter-filament distance')
    axes[2].set_ylabel('H - H_min')
    axes[2].set_title('Effective Potential Well')
    axes[2].grid(True, alpha=0.3)
    
    # Mark the equilibrium distance (minimum energy)
    eq_idx = np.argmin(ham_smooth)
    axes[2].axvline(x=distances[eq_idx], color='green', linestyle='--', 
                    alpha=0.7, label=f'Equilibrium d={distances[eq_idx]:.1f}')
    axes[2].legend()
    
    plt.tight_layout()
    plt.show()
    
    # Print summary
    print(f"Equilibrium distance: {distances[eq_idx]:.2f}")
    print(f"Minimum energy: {ham_smooth[eq_idx]:.2f}")
    print(f"Energy range: {np.min(ham_smooth):.2f} to {np.max(ham_smooth):.2f}")
    print(f"Max attractive force: {np.min(forces):.4f}")
    print(f"Max repulsive force: {np.max(forces):.4f}")
    
    return fig, (distances, mean_hams, std_hams, ham_smooth, forces)


# Example usage
if __name__ == "__main__":
    dirname = "20251210/prewet_enhancers_32_1-16"
    keya_values = list(range(1, 17))  # Adjust based on your actual keya range

    fig, results = plot_effective_force(
        dirname, 
        keya_values,
        keyb=1,
        burn_in_frac=0.5,      # Discard first 50% as burn-in
        n_samples=20,           # Take 20 samples from equilibrated region
        thinning=50,            # 50 frames between samples
        smoothing_window=3      # Smooth over 3 points
    )

    distances, mean_hams, std_hams, ham_smooth, forces = results
