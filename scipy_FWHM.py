import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
import os

def FindFWHM(x, y, peakIndex=0, height=0, distance=1, sortPeaks=True):
    """Find FWHM of a peak using scipy.signal functions"""
    
    # Find index of peaks
    peaks, _ = find_peaks(y, height=height, distance=distance)
    
    if len(peaks) == 0:
        raise ValueError("No peaks found in the data")
    # Sort the peaks by their height in descending order.
    if sortPeaks:
        peaks = peaks[np.argsort(y[peaks])][::-1]
    
    if peakIndex >= len(peaks):
        raise ValueError("Invalid peakIndex")
    
    peak = peaks[peakIndex]
    peak_value = y[peak]
    
    # Calculate width at half maximum
    widths, _, left_ips, right_ips = peak_widths(y, [peak], rel_height=0.5)
    
    # Get exact x-coordinates
    x_left = np.interp(left_ips[0], np.arange(len(x)), x)
    x_right = np.interp(right_ips[0], np.arange(len(x)), x)
    
    fwhm = x_right - x_left
    return fwhm, x_left, x_right, peak_value

def analyze_fwhm(file_path):
    """Calculate FWHM for each axis using the activity data"""
    # Load and reshape 3D data
    with h5py.File(file_path, 'r') as f:
        activities = f['activities'][()]
        voxel_config = eval(f.attrs['voxels'])
        nx, ny, nz = voxel_config['nXYZ']
        activities_3d = activities.reshape((nx, ny, nz))

    # Get center positions and axis arrays
    center_x, center_y, center_z = nx//2, ny//2, nz//2
    x_axis = np.linspace(-15, 15, nx)
    y_axis = np.linspace(-15, 15, ny)
    z_axis = np.linspace(-15, 15, nz)
    
    # Extract 1D activities through center
    x_activity = activities_3d[:, center_y, center_z]
    y_activity = activities_3d[center_x, :, center_z]
    z_activity = activities_3d[center_x, center_y, :]
    
    # Calculate FWHM for each axis
    try:
        fwhm_x, _, _, _ = FindFWHM(x_axis, x_activity)
        fwhm_y, _, _, _ = FindFWHM(y_axis, y_activity)
        fwhm_z, _, _, _ = FindFWHM(z_axis, z_activity)
    except ValueError as e:
        print(f"Error calculating FWHM for {file_path}: {e}")
        return np.nan, np.nan, np.nan

    return fwhm_x, fwhm_y, fwhm_z

"""
Process multiple epochs and create visualizations.
"""
def process_directories():
    base_path = '/fast_scratch_1/TIIGR/samuel/test'
    directories = ['ORIG_OUTPUT', 'OUTPUT_0_0_6', 'OUTPUT_0_0_0', 'OUTPUT_0_0_3', 'OUTPUT_0_0_4','OUTPUT_0_0_5', 'OUTPUT_0_0_7', 'OUTPUT_0_0_8','OUTPUT_kde_01_0_0_0','OUTPUT_kde_03_0_0_0']
    
    bandwidth_map = {
        'ORIG_OUTPUT': 'Original System Matrix',
        'OUTPUT_0_0_0': '0.030',
        'OUTPUT_0_0_3': '0.027',
        'OUTPUT_0_0_4': '0.033',
        'OUTPUT_0_0_5': '0.039',
        'OUTPUT_0_0_6': '0.021',
        'OUTPUT_0_0_7': '0.051',
        'OUTPUT_0_0_8': '0.066',
        'OUTPUT_kde_01_0_0_0': '0.1',
        'OUTPUT_kde_03_0_0_0': '0.3'
    }
    
    # Initialize results file
    csv_path = os.path.join(base_path, 'fwhm_results.csv')
    with open(csv_path, 'w') as f:
        f.write('Directory,Epoch,FWHM_X,FWHM_Y,FWHM_Z\n')
    
    for directory in directories:
        dir_path = os.path.join(base_path, directory)
        if not os.path.exists(dir_path):
            print(f"Skipping {directory} - directory not found")
            continue
            
        print(f"Processing {directory}...")
        epochs = range(1, 61)
        fwhm_x_values = []
        fwhm_y_values = []
        fwhm_z_values = []
        
        # Process each epoch
        for epoch in epochs:
            file_path = os.path.join(dir_path, f"{epoch:04d}.hdf5")
            if os.path.exists(file_path):
                fwhm_x, fwhm_y, fwhm_z = analyze_fwhm(file_path)
                fwhm_x_values.append(fwhm_x)
                fwhm_y_values.append(fwhm_y)
                fwhm_z_values.append(fwhm_z)
                
                # Save results to CSV
                with open(csv_path, 'a') as f:
                    f.write(f'{directory},{epoch},{fwhm_x:.3f},{fwhm_y:.3f},{fwhm_z:.3f}\n')
        
        plot_results(epochs, fwhm_x_values, fwhm_y_values, fwhm_z_values, directory, base_path)
    
    # Generate all comparison plots
    plot_fwhm_vs_bandwidth(directories, base_path, bandwidth_map)
    plot_yz_slices_comparison(directories, base_path, bandwidth_map)

"""
Create standardized visualization of FWHM evolution over epochs.
"""
def plot_results(epochs, fwhm_x, fwhm_y, fwhm_z, directory, base_path):
    """Create standardized visualization of FWHM evolution over epochs"""
    

    bandwidth_map = {
        'ORIG_OUTPUT': 'Original System Matrix', 
        'OUTPUT_0_0_2': '0.01',
        'OUTPUT_0_0_0': '0.030',
        'OUTPUT_0_0_3': '0.027',
        'OUTPUT_0_0_4': '0.033',
        'OUTPUT_0_0_5': '0.039',
        'OUTPUT_0_0_6': '0.021',
        'OUTPUT_0_0_7': '0.051',
        'OUTPUT_0_0_8': '0.066',
        'OUTPUT_kde_01_0_0_0': '0.1',
        'OUTPUT_kde_03_0_0_0': '0.3'
        
    }
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    
    # Set common y-axis limits and grid properties
    y_min = 0
    y_max = 14
    
    # Common plotting parameters
    plot_params = {
        'marker': 'o',
        'markersize': 4,
        'linestyle': 'None',  # Remove lines between points
        'linewidth': 1
    }
    
    # X-axis FWHM
    ax1.plot(epochs, fwhm_x, color='blue', **plot_params)
    ax1.set_title('X-axis')
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('FWHM')
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_ylim(y_min, y_max)
    ax1.set_xlim(0, 61)  
    
    # Y-axis FWHM
    ax2.plot(epochs, fwhm_y, color='green', **plot_params)
    ax2.set_title('Y-axis')
    ax2.set_xlabel('Epoch')
    ax2.set_ylabel('FWHM')
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.set_ylim(y_min, y_max)
    ax2.set_xlim(0, 61)  
    
    # Z-axis FWHM 
    ax3.plot(epochs, fwhm_z, color='cyan', **plot_params)
    ax3.set_title('Z-axis')
    ax3.set_xlabel('Epoch')
    ax3.set_ylabel('FWHM')
    ax3.grid(True, linestyle='--', alpha=0.7)
    ax3.set_ylim(y_min, y_max)
    ax3.set_xlim(0, 61)  
    
    # Add markers for last epoch
    last_epoch = epochs[-1]
    if directory == 'ORIG_OUTPUT':
        # Add crosses with "Double peak" label for ORIG_OUTPUT
        ax1.plot(last_epoch, fwhm_x[-1], 'rx', markersize=10, markeredgewidth=2, label='Multiple peaks')
        ax2.plot(last_epoch, fwhm_y[-1], 'rx', markersize=10, markeredgewidth=2, label='Multiple peaks')
        ax3.plot(last_epoch, fwhm_z[-1], 'rx', markersize=10, markeredgewidth=2, label='Multiple peaks')
        ax1.legend(loc='upper right')
    elif float(bandwidth_map[directory]) > 0.03:
        # Add crosses with "Oversmoothing" label for oversmoothed cases
        ax1.plot(last_epoch, fwhm_x[-1], 'rx', markersize=10, markeredgewidth=2, label='Oversmoothing')
        ax2.plot(last_epoch, fwhm_y[-1], 'rx', markersize=10, markeredgewidth=2, label='Oversmoothing')
        ax3.plot(last_epoch, fwhm_z[-1], 'rx', markersize=10, markeredgewidth=2, label='Oversmoothing')
        ax1.legend(loc='upper right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(base_path, f'fwhm_{directory}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def plot_fwhm_vs_bandwidth(directories, base_path, bandwidth_map):
    """Create plot of final FWHM values versus bandwidth with separate subplots for each axis"""
    # Skip ORIG_OUTPUT and collect data
    bandwidths = []
    final_fwhm_x = []
    final_fwhm_y = []
    final_fwhm_z = []
    
    for directory in directories:
        if directory == 'ORIG_OUTPUT':
            continue
            
        bandwidth = float(bandwidth_map[directory])
        
        # Get the last epoch's FWHM values
        file_path = os.path.join(base_path, directory, f"{60:04d}.hdf5")
        if os.path.exists(file_path):
            fwhm_x, fwhm_y, fwhm_z = analyze_fwhm(file_path)
            bandwidths.append(bandwidth)
            final_fwhm_x.append(fwhm_x)
            final_fwhm_y.append(fwhm_y)
            final_fwhm_z.append(fwhm_z)
    
    # Sort data points by bandwidth for proper curve plotting
    sorted_indices = np.argsort(bandwidths)
    bandwidths = np.array(bandwidths)[sorted_indices]
    final_fwhm_x = np.array(final_fwhm_x)[sorted_indices]
    final_fwhm_y = np.array(final_fwhm_y)[sorted_indices]
    final_fwhm_z = np.array(final_fwhm_z)[sorted_indices]
    
    # Create figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # Define specific bandwidth ticks (removed 0.01 and 0.046)
    bandwidth_ticks = [0.021, 0.027, 0.03, 0.033, 0.039, 0.051, 0.066, 0.1, 0.3]
    
    # Common plotting parameters
    optimal_bandwidth = 0.03
    y_min = min(min(final_fwhm_x), min(final_fwhm_y), min(final_fwhm_z)) - 0.5
    y_max = max(max(final_fwhm_x), max(final_fwhm_y), max(final_fwhm_z)) + 0.5
    
    # Common axis settings function
    def setup_axis(ax, data, color, title):
        ax.semilogx(bandwidths, data, f'{color}o', markersize=8)
        ax.set_title(f'{title} FWHM')
        ax.set_xlabel('Bandwidth')
        ax.set_ylabel('FWHM')
        # Add only horizontal grid lines
        ax.grid(True, which='both', linestyle='--', alpha=0.7, axis='y')
        ax.set_xscale('log')
        ax.set_xticks(bandwidth_ticks)
        ax.set_xticklabels(['{:.3f}'.format(x) for x in bandwidth_ticks], rotation=90)
        ax.set_ylim(y_min, y_max)
        # Add vertical line at optimal bandwidth
        ax.axvline(x=optimal_bandwidth, color='gray', linestyle='--', alpha=0.5, label='Optimal bandwidth')
    
    # Plot all three axes
    setup_axis(ax1, final_fwhm_x, 'b', 'X-axis')
    setup_axis(ax2, final_fwhm_y, 'g', 'Y-axis')
    setup_axis(ax3, final_fwhm_z, 'm', 'Z-axis')
    
    # Add red crosses for oversmoothing points
    for i, bw in enumerate(bandwidths):
        if bw > optimal_bandwidth and bw>0.039:
            ax1.plot(bw, final_fwhm_x[i], 'rx', markersize=10, markeredgewidth=2)
            ax2.plot(bw, final_fwhm_y[i], 'rx', markersize=10, markeredgewidth=2)
            ax3.plot(bw, final_fwhm_z[i], 'rx', markersize=10, markeredgewidth=2)
    
    # Add legend to first plot only
    handles = [
        plt.Line2D([0], [0], marker='o', color='b', linestyle='None', markersize=8, label='FWHM'),
        plt.Line2D([0], [0], color='gray', linestyle='--', label='Optimal bandwidth'),
        plt.Line2D([0], [0], marker='x', color='red', linestyle='None', markersize=10, 
                  markeredgewidth=2, label='Oversmoothing')
    ]
    ax1.legend(handles=handles, loc='upper right')
    
    # Adjust layout to prevent overlapping
    plt.subplots_adjust(bottom=0.2)  # Make room for rotated x-axis labels
    plt.tight_layout()
    
    plt.savefig(os.path.join(base_path, 'fwhm_vs_bandwidth.png'), dpi=300, bbox_inches='tight')
    plt.close()

def plot_yz_slices_comparison(directories, base_path, bandwidth_map):
    """Create YZ plane slices comparison for all bandwidths at the last epoch"""
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    # Skip ORIG_OUTPUT and collect data
    bandwidths = []
    yz_slices = []
    
    for directory in directories:
        if directory == 'ORIG_OUTPUT':
            continue
            
        bandwidth = float(bandwidth_map[directory])
        
        # Get the last epoch's data
        file_path = os.path.join(base_path, directory, "0060.hdf5")
        if os.path.exists(file_path):
            with h5py.File(file_path, 'r') as f:
                activities = f['activities'][()]
                voxel_config = eval(f.attrs['voxels'])
                nx, ny, nz = voxel_config['nXYZ']
                activities_3d = activities.reshape((nx, ny, nz))
                
                # Get YZ slice at X=0 (center)
                center_x = nx // 2  # This should be 24 for 48x48x48 grid
                yz_slice = activities_3d[center_x, :, :]
                
                bandwidths.append(bandwidth)
                yz_slices.append(yz_slice)
    
    # Sort by bandwidth
    sorted_indices = np.argsort(bandwidths)
    bandwidths = np.array(bandwidths)[sorted_indices]
    yz_slices = np.array(yz_slices)[sorted_indices]
    
    # Create figure with adjusted size and spacing
    nrows, ncols = 3, 3
    fig = plt.figure(figsize=(15, 12))  
    
    # Create gridspec with adjusted spacing
    gs = gridspec.GridSpec(nrows, ncols, figure=fig, 
                          hspace=0.4,  
                          wspace=0.6,  
                          left=0.1, right=0.9,  
                          bottom=0.1, top=0.9)
    
    # Plot each slice
    for idx, (bandwidth, yz_slice) in enumerate(zip(bandwidths, yz_slices)):
        if idx >= nrows * ncols:
            break
            
        row = idx // ncols
        col = idx % ncols
        ax = fig.add_subplot(gs[row, col])
        
        # Use individual color scale for each plot
        im = ax.imshow(yz_slice.T, origin='lower', cmap='viridis',
                      extent=[-15, 15, -15, 15])
        
        # Add colorbar for each subplot
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax, label='Activity')
        
        # Only show y-axis labels for leftmost plots
        if col == 0:
            ax.set_ylabel('Z position (mm)')
        else:
            ax.set_ylabel('')
            ax.set_yticklabels([])
        
        # Only show x-axis labels for bottom plots
        if row == nrows - 1:
            ax.set_xlabel('Y position (mm)')
        else:
            ax.set_xlabel('')
            ax.set_xticklabels([])
        
        ax.set_title(f'Bandwidth = {bandwidth:.3f}\n-1.0mm < X < 0.0mm', 
                    fontsize=10, pad=10)
    
    
    plt.savefig(os.path.join(base_path, 'yz_slices_comparison.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    process_directories()












