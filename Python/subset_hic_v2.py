import os
import numpy as np
import pandas as pd
import hicstraw
import csv
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec

def read_hic(hic_file):
    """
    Reads a HiC file using hicstraw and returns the HiCFile object.

    Parameters:
    - hic_file: Path to the .hic file.

    Returns:
    - hic: HiCFile object.
    """
    hic = hicstraw.HiCFile(hic_file)
    return hic

def get_contact(hic_file, regionX, regionY, resolution, output):
    """
    Fetches the contact data for specified regions from the .hic file and saves it to a CSV file.

    Parameters:
    - hic_file: Path to the .hic file.
    - regionX: Genomic region X in the format 'chr:start:end'.
    - regionY: Genomic region Y in the format 'chr:start:end'.
    - resolution: Resolution in base pairs.
    - output: Path to the output CSV file.

    Returns:
    - df: DataFrame containing the contact data.
    """
    # Fetch the counts and bin information to a list
    result = hicstraw.straw('observed', 'KR', hic_file, regionX, regionY, 'BP', resolution)
    
    # Store the retrieved data in a list of dictionaries
    data = []
    for entry in result:
        data.append({
            "region1": entry.binX,
            "region2": entry.binY,
            "counts": entry.counts
        })

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(data)

    # Save the DataFrame to a CSV file
    df.to_csv(output, index=False)

    print(f"Data saved to CSV file: {output}")
    return df

def get_matrix(hic_file, chr1, chr2, start1, end1, start2, end2, resolution):
    """
    Retrieves a subset of the HiC contact matrix for specified regions from the .hic file.

    Parameters:
    - hic_file: Path to the .hic file.
    - chr1: Chromosome 1.
    - chr2: Chromosome 2.
    - start1: Start position on chromosome 1.
    - end1: End position on chromosome 1.
    - start2: Start position on chromosome 2.
    - end2: End position on chromosome 2.
    - resolution: Resolution in base pairs.

    Returns:
    - numpy_matrix: Numpy array containing the contact matrix.
    """
    # Read the hic file
    hic = hicstraw.HiCFile(hic_file)

    # Subset the hic contact data
    mzd = hic.getMatrixZoomData(chr1, chr2, "observed", "KR", "BP", resolution)
    numpy_matrix = mzd.getRecordsAsMatrix(start1, end1, start2, end2)
    return numpy_matrix

def plot_heatmap(hic_file, chr1, chr2, start1, end1, start2, end2, resolution, output=None):
    """
    Plots a heatmap of the contact matrix for specified regions from the .hic file.

    Parameters:
    - hic_file: Path to the .hic file.
    - chr1: Chromosome 1.
    - chr2: Chromosome 2.
    - start1: Start position on chromosome 1.
    - end1: End position on chromosome 1.
    - start2: Start position on chromosome 2.
    - end2: End position on chromosome 2.
    - resolution: Resolution in base pairs.
    - output: Path to the output PNG file. If None, the plot will not be saved.

    Returns:
    - None
    """
    # Read the .hic file
    hic = hicstraw.HiCFile(hic_file)

    # Prepare the matrix
    matrix_object = hic.getMatrixZoomData(chr1, chr2, "observed", "KR", "BP", resolution)
    numpy_matrix = matrix_object.getRecordsAsMatrix(start1, end1, start2, end2)

    # Set the heatmap color
    REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])
    
    # Real plot function
    plt.matshow(numpy_matrix, cmap=REDMAP, vmin=0, vmax=30)

    # Save plot before plotting
    if output is not None:
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"Plot saved to PNG file: {output}")

    # Show the plot
    plt.show()

