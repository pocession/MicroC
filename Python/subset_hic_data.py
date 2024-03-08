import os
import numpy as np
import hicstraw
import math
import argparse

def parse_arguments():
    """
    Parse command line arguments.
    # python subset_hic_data.py --inputDir ./data/input.hiC ./Results/output.csv --outputDir --startChr chr2 --endChr chr2 
    # --startPosStartChr 112735986 --endPosStartChr 113204585 --startPosEndChr 112735986 --endPosEndChr 112735986 --res 5000
    """
    parser = argparse.ArgumentParser(description='Subset HiC Data.')
    parser.add_argument('--inputDir', type=str, required=True, help='the path of input hic file')
    parser.add_argument('--outputDir', type=str, required=True, help='the path of output file')
    parser.add_argument('--startChr', type=str, required=True, help='Start chromosome')
    parser.add_argument('--endChr', type=str, required=True, help='End chromosome')
    parser.add_argument('--startPosStartChr', type=int, required=True, help='Start position in start chromosome')
    parser.add_argument('--endPosStartChr', type=int, required=True, help='End position in start chromosome')
    parser.add_argument('--startPosEndChr', type=int, required=True, help='Start position in end chromosome')
    parser.add_argument('--endPosEndChr', type=int, required=True, help='End position in end chromosome')
    parser.add_argument('--res', type=int, required=True, help='Resolution')
    return parser.parse_args()

def subset_hic_data(inputDir, outputDir, startChr, endChr, startPosStartChr, endPosStartChr, startPosEndChr, endPosEndChr, res):
    """
    Subset HiC data based on the provided chromosome positions and resolution.
    """

    wd = os.getcwd()
    hicfile = wd + inputDir
    hic = hicstraw.HiCFile(hicfile)

    mzd = hic.getMatrixZoomData(startChr, endChr, "observed", "KR", "BP", res)
    numpy_matrix = mzd.getRecordsAsMatrix(startPosStartChr, endPosStartChr, startPosEndChr, endPosEndChr)

    output_dir = wd + outputDir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_name_with_extension = hicfile.split('/')[-1]
    outputfile = file_name_with_extension.rsplit('.', 1)[0]
    outputfile = output_dir + outputfile + "_extracted.csv"

    np.savetxt(outputfile, numpy_matrix, delimiter=",")
     # Print successful message
    print("Done! The file is save in " + output_dir)

def main():
    # Parse command line arguments
    args = parse_arguments()

    # Subset HiC data based on the provided arguments
    subset_hic_data(args.inputDir, args.outputDir, args.startChr, args.endChr, args.startPosStartChr, args.endPosStartChr, args.startPosEndChr, args.endPosEndChr, args.res)

if __name__ == "__main__":
    main()
