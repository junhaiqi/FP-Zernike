
# import zernike_version2
import sys
import numpy as np
import pickle

def extractInfoFromPklfile(pkl_file):
    
    f = open(pkl_file, 'rb')
    data = pickle.load(f)

    f.close()
    
    return data[list(data.keys())[0]]

def byInvFileGetDist(inv_file1, inv_file2):
    zernike_descriptor1 = extractInfoFromPklfile(inv_file1)
    zernike_descriptor2 = extractInfoFromPklfile(inv_file2)
    
    dist = np.linalg.norm(np.array(zernike_descriptor1)-np.array(zernike_descriptor2))

    return dist

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print('usage: python mainZernikeDist.py *.pkl *.pkl')
        print('Note: *.pkl only include one zernike discritior!')
        
    else:
        dist = byInvFileGetDist(inv_file1=sys.argv[1], inv_file2=sys.argv[2])
        print('The Zernike distance between %s and %s is: %f.'%(sys.argv[1], sys.argv[2], dist))
        print('Note: When the distance between descriptors is less than [PM:3, ATOM and PS: 2.5], it means that there is a high probability of strong similarity.')
