import pickle
import sys
import os
import numpy as np
import argparse

if 'module' not in sys.path:
    sys.path.append('module')
    
from calZernikeDescriptor import get_zernike
from scipy.spatial.distance import cdist

def extractInfoFromPklfile(pkl_file):
    """Get 3D Zernike descriptor from a 'pkl' file."""
    f = open(pkl_file, 'rb')
    data = pickle.load(f)
    f.close()
    return data

def getQueryZernike(pdbFile = '', mode = 'ATOM', prfixName = ''):
    # prfixName = pdbFile.split('/')[-1].split('.')[0].upper()
    zernikeData = {}
    if mode == 'ATOM':
        data = get_zernike(file_source = 'atom', input = pdbFile, gmmCount = 50, mode = 'atom', output = prfixName + '_ATOM.pkl')
        # zernikeData[prfixName] = data
    elif mode == 'PM':
        data = get_zernike(file_source = 'pymol', input = pdbFile, gmmCount = 50, mode = 'mesh', output = prfixName + '_PM.pkl')
        # zernikeData[prfixName] = data
    elif mode == 'PS':
        data = get_zernike(file_source = 'pymol', input = pdbFile, gmmCount = 50, mode = 'surface', output = prfixName + '_PS.pkl')
        # zernikeData[prfixName] = data
    elif mode == 'GM':
        data = get_zernike(file_source = 'gmconvert', input = pdbFile, gmmCount = 50, mode = 'surface', output = prfixName + '_GM.pkl') 
        # zernikeData[prfixName] = data
    else:
        print('Please input ATOM, PM, PS or GM!')
        sys.exit(-1)
    return data

def retrieveSystem(pdbFile = '', Zernikemode = 'ATOM', retrievalMode = 'SingleChain', outTopNum = 500, update = 'False', resFile = ''):
    prefixName = pdbFile.split('/')[-1].split('.')[0]
    
    if len(prefixName) == 5:
        prefixName = prefixName[0:4].upper() + prefixName[4] # pdb ID + chain ID, e.g., 1MBNA. Length = 5
    else:
        prefixName = prefixName.upper()

    dataFile = 'database/' + Zernikemode + 'All' + retrievalMode + 'Zernike.pkl'
    allZernikeData = extractInfoFromPklfile(dataFile)

    queryZernike = np.array([])
    if prefixName in allZernikeData.keys():
        queryZernike = allZernikeData[prefixName]
    else:
        if '.pdb' not in pdbFile:
            os.system('wget -P tempoutput/ https://files.rcsb.org/download/%s.pdb'%prefixName)
            os.system('mv tempoutput/%s.pdb tempoutput/%s.pdb'%(prefixName, prefixName.upper()))
            pdbFile = 'tempoutput/%s.pdb'%prefixName
        
        queryZernike = getQueryZernike(pdbFile = pdbFile, mode = Zernikemode, prfixName = prefixName)
        if update == 'True':
            os.system('rm -r %s'%dataFile)
            file = open(dataFile, 'ab')
            pickle.dump(allZernikeData, file, -1)
            file.close()

    keys = list(allZernikeData.keys())
    zernikes = np.array([allZernikeData[key] for key in allZernikeData])
    
    distList = (np.sum((((zernikes - queryZernike) ** 2)), axis = 1)) ** 0.5
    sortedList = sorted(distList)
    sortedDistListIndex = sorted(range(len(distList)), key=lambda k: distList[k], reverse=False)
    retrieveRes = [keys[item] for item in sortedDistListIndex][0:outTopNum]
    
    file = open(resFile, 'w')
    for i in range(len(retrieveRes)):
        file.write('%s %.2f\n'%(str(retrieveRes[i]), sortedList[i]))
    file.close()
    return retrieveRes

def initialization_parameters():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', type=str, required=True,
                        help='A path to the query protein sturcture(format: *.pdb).')
    parser.add_argument('-dm', type=str, required=False, default="ATOM",
                        help='The mode of zernike descriptor database. There are 3 modes, namely "ATOM", "PM" and "PS".')
    parser.add_argument('-otn', type=int, required=False, default=100,
                        help='The number of structures in the final output.')
    parser.add_argument('-op', type=str, required=True, \
                        help='An output file that records retrieve results.')
    parser.add_argument('-up', type=str, required=False, default="False",
                        help='A parameter that decides whether to update the database. \
                        If the parameter is set to "True" and the query structure is not in the database, \
                        a descriptor for this structure will be added to the database.')

    args = parser.parse_args()
    return args

def main():
    args = initialization_parameters()
    retrieveSystem(pdbFile = args.q, \
         Zernikemode = args.dm, \
         retrievalMode = 'SingleChain', \
         outTopNum = args.otn, \
         update = args.up, \
         resFile = args.op)

if __name__ == "__main__":
    main()
