import os
import numpy as np
import hicstraw
import math
import argparse

def parse_arguments():
    """
    Parse command line arguments.
    python ./Python/subset_hic_data.py --inputDir /data/44113_C_ctrl_43616_mc7contact_map.hic --outputDir /Results/processing/ --rowChr chr2 --columnChr chr2 --rowStart 112735986 --rowEnd 113204585 --columnStart 112735986 --columnEnd 113204585 --res 5000
    """
    parser = argparse.ArgumentParser(description='Subset HiC Data.')
    parser.add_argument('--inputDir', type=str, required=True, help='the path and name of input hic file')
    parser.add_argument('--outputDir', type=str, required=True, help='the path of output file')
    parser.add_argument('--rowChr', type=str, required=True, help='Start chromosome')
    parser.add_argument('--columnChr', type=str, required=True, help='End chromosome')
    parser.add_argument('--rowStart', type=int, required=True, help='Start position in start chromosome')
    parser.add_argument('--rowEnd', type=int, required=True, help='End position in start chromosome')
    parser.add_argument('--columnStart', type=int, required=True, help='Start position in end chromosome')
    parser.add_argument('--columnEnd', type=int, required=True, help='End position in end chromosome')
    parser.add_argument('--res', type=int, required=True, help='Resolution')
    return parser.parse_args()

def subset_hic_data(inputDir, outputDir, rowChr, columnChr, rowStart, rowEnd, columnStart, columnEnd, res):
    """
    Subset HiC data based on the provided chromosome positions and resolution.
    """

    wd = os.getcwd()
    hicfile = wd + inputDir
    hic = hicstraw.HiCFile(hicfile)

    mzd = hic.getMatrixZoomData(rowChr, columnChr, "observed", "KR", "BP", res)
    numpy_matrix = mzd.getRecordsAsMatrix(rowStart, rowEnd, columnStart, columnEnd)

    output_dir = wd + outputDir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_name_with_extension = hicfile.split('/')[-1]
    outputfile = file_name_with_extension.rsplit('.', 1)[0]
    outputfile = output_dir + outputfile + "_extracted.csv"

    np.savetxt(outputfile, numpy_matrix, delimiter=",")
     # Print successful message
    print("Done! The interaction matrix file is saved in " + outputfile)

    return(numpy_matrix)
def generate_startPos_index(numpy_matrix, input_dir, output_dir, rowStart, rowEnd, columnStart, columnEnd, res):
    """
    Generate the position index for the interaction matrix and write into a txt file
    The position index of each tile is the midpoint of that tile
    For example, 5000 stands for the tile of 2500-7500 when resolution = 5000
    The length of the last tile may not be divisible by resolution
    In this case, the interval between the last pos index and last-1 pos index would be smaller than resolution
    """

    ## Get the matrix dimension
    rowNum = numpy_matrix.shape[0]
    columnNum = numpy_matrix.shape[1]

    ## Create the rowIndex
    rowIndex = []
    newRowStart = rowStart+math.ceil(res/2)
    rowIndex.append(newRowStart)
    for i in range(1,(rowNum)):
        newRowStart = rowIndex[i-1] + res
        rowIndex.append(newRowStart)

    rowLastInterval = math.ceil((rowEnd - rowIndex[(len(rowIndex)-2)])/2)
    rowIndex[len(rowIndex)-1] = rowIndex[len(rowIndex)-2] + rowLastInterval

    ## Create the columnIndex
    columnIndex = []
    newColumnStart = columnStart+math.ceil(res/2)
    columnIndex.append(newColumnStart)
    for i in range(1,(columnNum)):
        newColumnStart = columnIndex[i-1] + res
        columnIndex.append(newColumnStart)

    columnLastInterval = math.ceil((columnEnd - columnIndex[(len(columnIndex)-2)])/2)
    columnIndex[len(columnIndex)-1] = columnIndex[len(columnIndex)-2] + columnLastInterval

    ## Generate the starting position index 
    ## row[i]-column[j] index
    startingPosIndex = []
    for i in range(len(rowIndex)):
        for j in range(len(columnIndex)):
            index = str(rowIndex[i]) + "-" + str(columnIndex[j])
            startingPosIndex.append(index)

    ## Get the output file name
    wd = os.getcwd()
    hicfile = wd + input_dir

    ## Extract the last part after the last slash, then remove the extension
    file_name_with_extension = hicfile.split('/')[-1]
    outputfile = file_name_with_extension.rsplit('.', 1)[0]
    outputfile = wd + output_dir + outputfile + "_PosIndex.txt"
    data = startingPosIndex
    with open(outputfile, "w") as txt_file:
        for line in data:
            txt_file.write(str(line) + "\n")
    
    print("Done! The position index file is saved in " + outputfile)

def main():
    # Parse command line arguments
    args = parse_arguments()

    # Subset HiC data based on the provided arguments
    numpy_matrix = subset_hic_data(args.inputDir, args.outputDir, args.rowChr, args.columnChr, args.rowStart, args.rowEnd, args.columnStart, args.columnEnd, args.res)

    # Generate the starting position index
    generate_startPos_index(numpy_matrix, args.inputDir, args.outputDir, args.rowStart, args.rowEnd, args.columnStart, args.columnEnd, args.res)

if __name__ == "__main__":
    main()
