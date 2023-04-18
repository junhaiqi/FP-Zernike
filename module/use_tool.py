
import os

current_file_path = os.path.dirname(os.path.abspath('tool'))
prefix = current_file_path + '/tool/'  # absolute path.


def get_pdbSurface(pdb_file, GMMCount=100):

    pre_name = pdb_file.split('/')[-1]
    pre_name = pre_name.split('.')[0]
    print('Generating the surface of %s' % pdb_file)
    os.system('%sgmconvert A2G -ipdb %s -ng %d -maxsize 64 -ogmm tempoutput/%s_ng%d.gmm >/dev/null' %
              (prefix, pdb_file, GMMCount, pre_name, GMMCount))
    os.system('%sgmconvert G2S -igmm tempoutput/%s_ng%d.gmm -owrl tempoutput/%s_ng%d_surface.wrl >/dev/null' %
              (prefix, pre_name, GMMCount, pre_name, GMMCount))
    print('End of generate!')
    
def useRNAAlignGetTMscoreRMSD(pdb_file1, pdb_file2):
    
    outDict = {}
    pre_name1 = pdb_file1.split('/')[-1]
    pre_name1 = pre_name1.split('.')[0]
    
    pre_name2 = pdb_file2.split('/')[-1]
    pre_name2 = pre_name2.split('.')[0]
    
    os.system('%sRNAalign %s %s -outfmt 2 >tempoutput/%s_%s_RNAalignOut.txt'
              %(prefix, pdb_file1, pdb_file2, pre_name1, pre_name2))
    
    with open('tempoutput/%s_%s_RNAalignOut.txt'%(pre_name1, pre_name2)) as f:
        lines = f.readlines()
        outLine = lines[1].strip('\n')
        
        outList = outLine.split('\t')
        
        resList = outList[2:5]
        
        if len(resList) == 3:
        
            outDict['TMScore'] = (float(resList[0]) + float(resList[1])) / 2
            outDict['RMSD'] = float(resList[2])
            return outDict
        
        else:
            outDict['TMScore'] = None
            outDict['RMSD'] = None
            
            return outDict
        # print(lines)

def useTMscore(pdbFile1 = '', pdbFile2 = ''):
    
    pre_name1 = pdbFile1.split('/')[-1]
    pre_name1 = pre_name1.split('.')[0]
    
    pre_name2 = pdbFile2.split('/')[-1]
    pre_name2 = pre_name2.split('.')[0]
    
    os.system('%sTMscore %s %s >tempoutput/%s_%s_TMscoreOut.txt'
              %(prefix, pdbFile1, pdbFile2, pre_name1, pre_name2))
    
    file_name = 'tempoutput/%s_%s_TMscoreOut.txt'%(pre_name1, pre_name2)
    fail_index = 'no'
    
    with open(file_name,'r') as c:
        lines = c.readlines()
        if lines == []:
            return fail_index
        
        for line in lines:
            if 'TM-score' in line and '=' in line:
                if 'Superposition' not in line:
                    line = line.strip('\n')
                    line = line.split('=')
                    line = line[1].split('(')
                    line[0] = float(line[0])
                    return line[0]
                
def useDTW(valueFile1 = '', valueFile2 = ''):
    
    pre_name1 = valueFile1.split('/')[-1]
    pre_name1 = pre_name1.split('.')[0]
    
    pre_name2 = valueFile2.split('/')[-1]
    pre_name2 = pre_name2.split('.')[0]
    
    os.system('%sdtw_zscore3 %s %s >tempoutput/%s_%s_DTWdist.txt'
              %(prefix, valueFile1, valueFile2, pre_name1, pre_name2))
    

    file_name = 'tempoutput/%s_%s_DTWdist.txt'%(pre_name1, pre_name2)
    
    dtwDist = [float(line) for line in open(file_name)][0]
    
    return dtwDist


def useDeepScore(pdbFile1 = '', pdbFile2 = ''):
    pre_name1 = pdbFile1.split('/')[-1]
    pre_name1 = pre_name1.split('.')[0]
    
    pre_name2 = pdbFile2.split('/')[-1]
    pre_name2 = pre_name2.split('.')[0]
    
    os.system('%sDeepScore %s %s >tempoutput/%s_%s_DeepscoreOut.txt'
              %(prefix, pdbFile1, pdbFile2, pre_name1, pre_name2))
    
    file_name = 'tempoutput/%s_%s_DeepscoreOut.txt'%(pre_name1, pre_name2)
    fail_index = 'no'
    
    with open(file_name,'r') as c:
        lines = c.readlines()
        if lines == []:
            return fail_index
        
        scoreIndex = 0
        for i in range(0, len(lines)):
            if '# BLOSUM CLESUM DeepScore SeqID LALI RMSD(A) TMscore MAXSUB GDT_TS GDT_HA uGDT' in lines[i]:
                scoreIndex = i + 1
                break
        
        scoreList = [item for item in lines[scoreIndex].split(' ') if item != '']
        
        RMSDScore = float(scoreList[5])
        TMscore = float(scoreList[6])
        
        return {'TMscore':TMscore, 'RMSD':RMSDScore}
    

def useGmfit(GmmFile1, GmmFile2):
    
    # prefixName1 = GmmFile1.split('/')[-1].split('_')[0]
    # prefixName2 = GmmFile2.split('/')[-1].split('_')[0]

    prefixName1 = GmmFile1.split('/')[-1][0:len(GmmFile1.split('/')[-1]) - 9]
    prefixName2 = GmmFile2.split('/')[-1][0:len(GmmFile2.split('/')[-1]) - 9]
    
    command = '%sgmfit -cg %s -sg1 %s -I P -ochimera T >tempoutput/%s_%s_gmfitscoreOut.txt' \
                %(prefix, GmmFile1, GmmFile2, prefixName1, prefixName2)
                
    os.system(command=command)
        
    
if __name__ == "__main__":
    test = useGmfit(GmmFile1='Benchmark/gmfit/SEQGMM/3cptA_g50.gmm', 
                    GmmFile2='Benchmark/gmfit/SEQGMM/3cptA_g50.gmm')

