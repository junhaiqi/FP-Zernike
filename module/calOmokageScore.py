import os
from statistics import mean
import numpy as np
import sys
import math
from sklearn.decomposition import PCA
import argparse

def initialization_parameters():
    parser = argparse.ArgumentParser()

    parser.add_argument('-pdb1', type=str, required=True,
                        help='input file, it is a pdb file.')
    parser.add_argument('-pdb2', type=str, required=True,
                        help='input file, it is a pdb file.')
    parser.add_argument('-info', type=str, required=False,
                        help='This is a script for calculating ReOmokage score between two structures.')

    args = parser.parse_args()

    return args

def list_remove_empty(list_temp):
    for i in list_temp[:]:
        if i == '':
            list_temp.remove(i)
    return list_temp

def L_2distance(list_a,list_b):
    if len(list_a) != len(list_b):
        print('List length error!')
        exit(-1)
    else:
        distance_2 = 0
        for i in range(0,len(list_a)):
            temp = list_a[i]-list_b[i]
            distance_2 = math.pow(temp,2) + distance_2     
    distance_sqrt = math.sqrt(distance_2)
    return distance_sqrt
            
def Get_list_distanceDistrbution(list_temp):
    distance_list = []
    for i in range(0,len(list_temp)-1):
        for j in range(i+1,len(list_temp)):
            distance = L_2distance(list_temp[i],list_temp[j])
            distance_list.append(distance)
    distance_list.sort()
    return distance_list

def PCA_simple(feature_list):
    pca = PCA(n_components=3)
    pca.fit(np.array(feature_list))
    temp_list = pca.fit_transform(np.array(feature_list))
    list_1 = []
    list_2 = []
    list_3 = []
    PCA_list = []
    
    for i in range(0,len(temp_list)):
        list_1.append(temp_list[i][0])
        list_2.append(temp_list[i][1])
        list_3.append(temp_list[i][2])
    list_1 = np.array(list_1)
    list_2 = np.array(list_2)
    list_3 = np.array(list_3)
    PCA_list.append(np.std(list_1))
    PCA_list.append(np.std(list_2))
    PCA_list.append(np.std(list_3))
    return PCA_list

def Out_25(feature_list):
    center = []
    mean_1 = 0
    mean_2 = 0
    mean_3 = 0
    for i in range(0,len(feature_list)):
        mean_1 = mean_1 + feature_list[i][0]
        mean_2 = mean_2 + feature_list[i][1]
        mean_3 = mean_3 + feature_list[i][2]
    l = len(feature_list)
    center.append(mean_1/l)
    center.append(mean_2/l)
    center.append(mean_3/l)
    distance_list = []
    out_list = []
    for item in feature_list:
        distance_list.append(L_2distance(item,center))
    # print(distance_list)
    sortd_id = sorted(range(len(distance_list)), key=lambda k : distance_list[k], reverse=True)
    for item in sortd_id:
        out_list.append(feature_list[item])
    return out_list[0:25]

def Get_feature_points(pdb_file, mode=50):
    # os.system('./quanpdb %s %sfeature_file.pdb'%(pdb_file,pdb_file[0:len(pdb_file)-4]))
    pdbID = pdb_file.split('/')[-1].split('.')[0]
    if mode == 50:
        try:
            temp_file = open('ReOmoFea/%sfeature_file50.pdb'%pdbID,'r')
        except:
            print('####################################################')
            print('Extracting 50 feature points of structure %s...'%pdb_file)
            os.system('./tool/extractFeaturePoints50 %s ReOmoFea/%sfeature_file50.pdb >/dev/null'%(pdb_file, pdbID))
            print('End of extraction!')
            print('####################################################\n')
            temp_file = open('ReOmoFea/%sfeature_file50.pdb'%pdbID,'r')
    elif mode == 30:
        try:
            temp_file = open('ReOmoFea/%sfeature_file30.pdb'%pdbID,'r')
        except:
            print('####################################################')
            print('Extracting 30 feature points of structure %s...'%pdb_file)
            os.system('./tool/extractFeaturePoints30 %s ReOmoFea/%sfeature_file30.pdb >/dev/null'%(pdb_file, pdbID))
            print('End of extraction!')
            print('####################################################\n')
            temp_file = open('ReOmoFea/%sfeature_file30.pdb'%pdbID,'r')

    feature_list = []
    for line in temp_file:
        line = line.strip('\n')
        line = line.split(' ')
        line = list_remove_empty(line)
        line = line[5:8] 
        new_line = []
        for i in range(0,len(line)):
            try:
                line[i] = float(line[i])
                new_line.append(line[i])
            except ValueError:
                test_list = []
                for j in range(0,len(line[i])):
                    if line[i][j] == '-':
                        test_list.append(j)
                if len(test_list) == 2: 
                    item_1 = line[i][test_list[0]:test_list[1]]
                    item_2 = line[i][test_list[1]:len(line[i])]
                    if i == 0:
                       new_line.append(float(item_1))
                       new_line.append(float(item_2))
                    if i == 1:
                        new_line.append(float(item_1))
                        new_line.append(float(item_2))
                if len(test_list) == 1:
                    item_1 = line[i][0:test_list[0]]
                    item_2 = line[i][test_list[0]:len(line[i])]
                    if i == 0:
                      new_line.append(float(item_1))
                      new_line.append(float(item_2))
                    if i == 1:
                        new_line.append(float(item_1))
                        new_line.append(float(item_2))
        feature_list.append(new_line[0:3])
    return feature_list
        
def GetDRprofile(size,feature_list,pdb_file):
    if size == 50:
        # pdb_file = sys.argv[1]
        DR_profile = Get_list_distanceDistrbution(feature_list)
        return DR_profile,feature_list
    if size == 30:
        # pdb_file = sys.argv[1]
        feature_list =  Get_feature_points(pdb_file, mode = 30)
        DR_profile = Get_list_distanceDistrbution(feature_list)
        return DR_profile,feature_list
    if size == 25:
        temp_list = Out_25(feature_list)
        DR_profile = Get_list_distanceDistrbution(temp_list)
        return DR_profile,temp_list
    if size == 'pca':
        PCA_list = PCA_simple(feature_list)
        return PCA_list

def ByDRprofileGetiDRprofile(DR_profile,win_len):
    iDR_profile = []
    L = len(DR_profile)
    for i in range(0,len(DR_profile)):
        temp = min(win_len,L-i)
        p_i = 0
        for t in range(0,temp):
            p_i += (DR_profile[t+i]-DR_profile[i])
        iDR_profile.append(p_i)
    return iDR_profile

def averageList(length = 24, list = [1]*1225):
    
    resList = [mean(list[i*length:(i+1)*length]) for i in range(0, int(len(list) / length))]
    
    if int(len(list) / length) * length < len(list):
    
        addValue = list[int(len(list) / length)*length:][0]
        
        resList.append(addValue)
    
    return resList

def SK_prepare(iDR_profile1,iDR_profile2):
    a = np.array(iDR_profile1)-np.array(iDR_profile2)
    b = np.array(iDR_profile1)+np.array(iDR_profile2)
    son = 0
    for item in a:
        son = son + math.pow(item,2)
    mother = 0
    for item in b:
        mother = mother + math.pow(item,2)
    sk = math.sqrt(son/mother)
    return 1-sk

def SK(K,pdb_file1,pdb_file2,feature_list_1,feature_list_2):
    if K == 30:
        DR_profile1 = GetDRprofile(30,feature_list_1,pdb_file1)[0]
        iDR_profile1 = ByDRprofileGetiDRprofile(DR_profile1,90) 
        iDR_profile1 = averageList(length = 4, list = iDR_profile1)
        
        DR_profile2 = GetDRprofile(30,feature_list_2,pdb_file2)[0]
        iDR_profile2 = ByDRprofileGetiDRprofile(DR_profile2,90)
        iDR_profile2 = averageList(length = 4, list = iDR_profile2)
        
        sk = SK_prepare(iDR_profile1,iDR_profile2)
        return sk
    if K == 50:
        DR_profile1 = GetDRprofile(50,feature_list_1,pdb_file1)[0]
        iDR_profile1 = ByDRprofileGetiDRprofile(DR_profile1,368)
        iDR_profile1 = averageList(length = 24, list = iDR_profile1)
        
        DR_profile2 = GetDRprofile(50,feature_list_2,pdb_file2)[0]
        iDR_profile2 = ByDRprofileGetiDRprofile(DR_profile2,368)
        iDR_profile2 = averageList(length = 24, list = iDR_profile2)
        
        sk = SK_prepare(iDR_profile1,iDR_profile2)
        return sk
    if K == 25:
        DR_profile1 = GetDRprofile(25,feature_list_1,pdb_file1)[0]
        iDR_profile1 = ByDRprofileGetiDRprofile(DR_profile1,131)
        iDR_profile1 = averageList(length = 4, list = iDR_profile1)
        
        DR_profile2 = GetDRprofile(25,feature_list_2,pdb_file2)[0]
        iDR_profile2 = ByDRprofileGetiDRprofile(DR_profile2,131)
        iDR_profile2 = averageList(length = 4, list = iDR_profile2)
        
        sk = SK_prepare(iDR_profile1,iDR_profile2)
        return sk
    if K == 'pca':
        PCA_list1 = GetDRprofile('pca',feature_list_1,pdb_file1)
        PCA_list2 = GetDRprofile('pca',feature_list_2,pdb_file2)
        sk = SK_prepare(PCA_list1,PCA_list2)
        return sk
    
def SRMS(pdb_file1,pdb_file2,feature_list_1,feature_list_2):
    a = SK(30,pdb_file1,pdb_file2,feature_list_1,feature_list_2)
    b = SK(50,pdb_file1,pdb_file2,feature_list_1,feature_list_2)
    c = SK(25,pdb_file1,pdb_file2,feature_list_1,feature_list_2)
    d = SK('pca',pdb_file1,pdb_file2,feature_list_1,feature_list_2)
    s = math.pow(a,2) + math.pow(b,2) + math.pow(c,2) + math.pow(d,2)
    s = math.sqrt(s/4)
    p = 2.17938653227261
    s = 2*math.pow(s,p)-1
    return s

def main():
    args = initialization_parameters()
    pdb_file1 = args.pdb1
    pdb_file2 = args.pdb2

    preName1 = pdb_file1.split('/')[-1].split('.')[0]
    preName2 = pdb_file2.split('/')[-1].split('.')[0]
    
    feature_list_1 = Get_feature_points(pdb_file1, mode = 50)
    feature_list_2 = Get_feature_points(pdb_file2, mode = 50)
    score = SRMS(pdb_file1,pdb_file2,feature_list_1,feature_list_2)
    print('ReOmokage score between %s and %s is:%f.'%(pdb_file1, pdb_file2, score))
    print('Note: When the score between structures is more than 0.8, it means that there is a high probability of strong similarity.')

    return ((preName1, preName2, score))

if __name__ == "__main__":
    main()