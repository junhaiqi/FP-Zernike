import sys
import pickle
import os

def extractInfoFromPklfile(pkl_file):
    """Get 3D Zernike descriptor from a 'pkl' file."""
    f = open(pkl_file, 'rb')
    data = pickle.load(f)
    f.close()
    keys = list(data.keys())
    return data[keys[0]]

def main(zernikeFloderPath = '', outFile = ''):

    resDict = {}
    pklList = os.listdir(path=zernikeFloderPath)
    for item in pklList:
        prefixName = item.split('.')[0].upper()
        data = extractInfoFromPklfile(zernikeFloderPath + '/' + item)
        resDict[prefixName] = data

    file = open(outFile, 'ab')
    pickle.dump(resDict, file, -1)
    file.close()


if __name__ == '__main__':
    # main(zernikeFloderPath = 'exampleAtomZernikeDescriptor', outFile = 'ATOMAllSingleChainZernike.pkl')

    main(zernikeFloderPath = 'examplePymolMeshZernikeDescriptor', outFile = 'PMAllSingleChainZernike.pkl')
    main(zernikeFloderPath = 'examplePymolSurfaceZernikeDescriptor', outFile = 'PSAllSingleChainZernike.pkl')
    main(zernikeFloderPath = 'examplegmconvertZernikeDescriptor', outFile = 'GMAllSingleChainZernike.pkl')