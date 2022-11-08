import os
import sys
from time import time
import argparse
from multiprocessing import Pool

def getPymolSurfaceZernike(pdbFile, pklDirPrefixName = 'Random100'):

    prefixName = pdbFile.split('/')[-1].split('.')[0]

    command = 'python mainZernikeCalculation.py -i %s -o \
    %sPymolSurfaceZernikeDescriptor/%s.pkl -file_source pymol -mode surface >>log/%sPymolSurface.txt' \
            %(pdbFile, pklDirPrefixName, prefixName, prefixName)

    os.system(command = command)

    os.system('rm tempoutput/%s_surface.wrl'%prefixName) # remove representation file for saving space.

def getPymolMeshZernike(pdbFile, pklDirPrefixName = 'Random100'):

    prefixName = pdbFile.split('/')[-1].split('.')[0]

    command = 'python mainZernikeCalculation.py -i %s -o \
    %sPymolMeshZernikeDescriptor/%s.pkl -file_source pymol -mode mesh >>log/%sPymolMesh.txt' \
            %(pdbFile, pklDirPrefixName, prefixName, prefixName)

    os.system(command = command)

    os.system('rm tempoutput/%s_mesh.wrl'%prefixName)

def getatomZernike(pdbFile, pklDirPrefixName = 'Random100'):

    prefixName = pdbFile.split('/')[-1].split('.')[0]
    command = 'python mainZernikeCalculation.py -i %s -o \
    %sAtomZernikeDescriptor/%s.pkl -file_source atom -mode atom >>log/%sAtom.txt' \
            %(pdbFile, pklDirPrefixName, prefixName, prefixName)

    os.system(command = command)

def getgmconvertZernike(pdbFile, pklDirPrefixName = 'Random100'):

    prefixName = pdbFile.split('/')[-1].split('.')[0]

    command = 'python mainZernikeCalculation.py -i %s -o \
    %sgmconvertZernikeDescriptor/%s.pkl -file_source gmconvert -mode surface >>log/%sgmconvert.txt' \
            %(pdbFile, pklDirPrefixName, prefixName, prefixName)

    os.system(command = command)

    os.system('rm tempoutput/%s_ng%d*'%(prefixName, 50))


def getZernikeDataBase(sourcePDBDir = 'Random100', threadNum = 4, mode=None):

    if sourcePDBDir[-1] == '/':
        sourcePDBDir = sourcePDBDir[0:-1]
        print(sourcePDBDir)
    if mode == None:
    
        os.system('mkdir %sgmconvertZernikeDescriptor'%sourcePDBDir)
        os.system('mkdir %sAtomZernikeDescriptor'%sourcePDBDir)
        os.system('mkdir %sPymolMeshZernikeDescriptor'%sourcePDBDir)
        os.system('mkdir %sPymolSurfaceZernikeDescriptor'%sourcePDBDir)

        start_time = time()
        successList = [item.split('.')[0] for item in os.listdir('%sPymolSurfaceZernikeDescriptor'%sourcePDBDir)]
        pdbList = ['%s/'%sourcePDBDir + item for item in os.listdir(sourcePDBDir) \
                   if '.pdb' in item and item.split('.')[0] not in successList]
        prefixNameList = [sourcePDBDir]*len(pdbList)
        args = [(pdbList[i], prefixNameList[i]) for i in range(len(pdbList))]
        pool = Pool(threadNum)
        pool.starmap(getPymolSurfaceZernike, args)
        pool.close()
        pool.join()
        end_time = time()
        # print('Total time to generate %sPymolSurfaceZernikeDescriptor: %fs'%(sourcePDBDir, (end_time - start_time)))

        start_time = time()
        successList = [item.split('.')[0] for item in os.listdir('%sPymolMeshZernikeDescriptor'%sourcePDBDir)]
        pdbList = ['%s/'%sourcePDBDir + item for item in os.listdir(sourcePDBDir) \
                   if '.pdb' in item and item.split('.')[0] not in successList]
        prefixNameList = [sourcePDBDir]*len(pdbList)
        args = [(pdbList[i], prefixNameList[i]) for i in range(len(pdbList))]
        pool = Pool(threadNum)
        pool.starmap(getPymolMeshZernike, args)
        pool.close()
        pool.join()
        end_time = time()
        # print('Total time to generate %sPymolMeshZernikeDescriptor: %fs'%(sourcePDBDir, (end_time - start_time)))

        start_time = time()
        successList = [item.split('.')[0] for item in os.listdir('%sgmconvertZernikeDescriptor'%sourcePDBDir)]
        pdbList = ['%s/'%sourcePDBDir + item for item in os.listdir(sourcePDBDir) \
                   if '.pdb' in item and item.split('.')[0] not in successList]
        prefixNameList = [sourcePDBDir]*len(pdbList)
        args = [(pdbList[i], prefixNameList[i]) for i in range(len(pdbList))]
        pool = Pool(threadNum)
        pool.starmap(getgmconvertZernike, args)
        pool.close()
        pool.join()
        end_time = time()
        # print('Total time to generate %sgmconvertZernikeDescriptor: %fs'%(sourcePDBDir, (end_time - start_time)))

        start_time = time()
        successList = [item.split('.')[0] for item in os.listdir('%sAtomZernikeDescriptor'%sourcePDBDir)]
        pdbList = ['%s/'%sourcePDBDir + item for item in os.listdir(sourcePDBDir) \
                   if '.pdb' in item and item.split('.')[0] not in successList]
        prefixNameList = [sourcePDBDir]*len(pdbList)
        args = [(pdbList[i], prefixNameList[i]) for i in range(len(pdbList))]
        pool = Pool(threadNum)
        pool.starmap(getatomZernike, args)
        pool.close()
        pool.join()
        end_time = time()
        # print('Total time to generate %sAtomZernikeDescriptor: %fs'%(sourcePDBDir, (end_time - start_time)))

    elif mode == 'GM':
        os.system('mkdir %sgmconvertZernikeDescriptor'%sourcePDBDir)

        successList = [item.split('.')[0] for item in os.listdir('%sgmconvertZernikeDescriptor'%sourcePDBDir)]
        pdbList = ['%s/'%sourcePDBDir + item for item in os.listdir(sourcePDBDir) \
                   if '.pdb' in item and item.split('.')[0] not in successList]
        prefixNameList = [sourcePDBDir]*len(pdbList)
        args = [(pdbList[i], prefixNameList[i]) for i in range(len(pdbList))]

        start_time = time()
        pool = Pool(threadNum)
        pool.starmap(getgmconvertZernike, args)
        pool.close()
        pool.join()
        end_time = time()
        # print('Total time to generate %sgmconvertZernikeDescriptor: %fs'%(sourcePDBDir, (end_time - start_time)))

    elif mode == 'PM':
        os.system('mkdir %sPymolMeshZernikeDescriptor'%sourcePDBDir)
        
        successList = [item.split('.')[0] for item in os.listdir('%sPymolMeshZernikeDescriptor'%sourcePDBDir)]
        pdbList = ['%s/'%sourcePDBDir + item for item in os.listdir(sourcePDBDir) \
                   if '.pdb' in item and item.split('.')[0] not in successList]
        prefixNameList = [sourcePDBDir]*len(pdbList)
        args = [(pdbList[i], prefixNameList[i]) for i in range(len(pdbList))]

        start_time = time()
        pool = Pool(threadNum)
        pool.starmap(getPymolMeshZernike, args)
        pool.close()
        pool.join()
        end_time = time()
        # print('Total time to generate %sPymolMeshZernikeDescriptor: %fs'%(sourcePDBDir, (end_time - start_time)))
    
    elif mode == 'PS':
        os.system('mkdir %sPymolSurfaceZernikeDescriptor'%sourcePDBDir)
        
        successList = [item.split('.')[0] for item in os.listdir('%sPymolSurfaceZernikeDescriptor'%sourcePDBDir)]
        pdbList = ['%s/'%sourcePDBDir + item for item in os.listdir(sourcePDBDir) \
                   if '.pdb' in item and item.split('.')[0] not in successList]
        prefixNameList = [sourcePDBDir]*len(pdbList)
        args = [(pdbList[i], prefixNameList[i]) for i in range(len(pdbList))]

        start_time = time()
        pool = Pool(threadNum)
        pool.starmap(getPymolSurfaceZernike, args)
        pool.close()
        pool.join()
        end_time = time()
        # print('Total time to generate %sPymolSurfaceZernikeDescriptor: %fs'%(sourcePDBDir, (end_time - start_time)))
    
    elif mode == 'ATOM':
        os.system('mkdir %sAtomZernikeDescriptor'%sourcePDBDir)
        
        successList = [item.split('.')[0] for item in os.listdir('%sAtomZernikeDescriptor'%sourcePDBDir)]
        pdbList = ['%s/'%sourcePDBDir + item for item in os.listdir(sourcePDBDir) \
                   if '.pdb' in item and item.split('.')[0] not in successList]
        prefixNameList = [sourcePDBDir]*len(pdbList)
        args = [(pdbList[i], prefixNameList[i]) for i in range(len(pdbList))]

        start_time = time()
        pool = Pool(threadNum)
        pool.starmap(getatomZernike, args)
        pool.close()
        pool.join()
        end_time = time()
        # print('Total time to generate %sAtomZernikeDescriptor: %fs'%(sourcePDBDir, (end_time - start_time)))

def initialization_parameters():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True,
                        help='A path to the folder that contains the pdb file.')
    parser.add_argument('-t', type=int, required=False, default=8,
                        help='Number of processes, default=8.')
    parser.add_argument('-mode', type=str, required=False, default=None,
                        help='Protein representation type. The following options are available: PM, PS, GM, ATOM, defalut=None.')
    parser.add_argument('-info', type=str, required=False, default=None,
                        help='This is a script to build the descriptor database.')

    args = parser.parse_args()
    return args

def main():
    args = initialization_parameters()

    getZernikeDataBase(sourcePDBDir = args.i, threadNum = args.t, mode = args.mode)


if __name__ == "__main__":
    main()
