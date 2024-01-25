# FP-Zernike: An open-source structural database construction toolkit for fast and accurate structure retrieval

The main functions of FP-Zernike as follows:

1. Given a protein(or RNA) structure(format: *.pdb), FP-Zernike can calculate its Zernike descriptors[4 types, namely PM-Zernike, PS-Zernike, ATOM-Zernike and GM-Zernike]. Based on the results of the experiment, we recommend the user to calculate PM-Zernike, PS-Zernike or ATOM-Zernike.

2. Given a customized dataset[only include *.pdb file], FP-Zernike can calculate all Zernike descriptors about it.

3. Given a query structure[single-chain], structure retrieval is done on a dataset containing 590685 protein chains[run time: ~4s-9s].

4. Based on the Euclidean distance of descriptors, the similarity of structures is preliminarily evaluated. When the Euclidean distance between descriptors is less than 2, it means that there is strong similarity.

5. Based on 'Romokage', the similarity between structures is measured.

FP-Zernike is an open-source toolkit, whose source code and related data are accessible at https://github.com/junhaiqi/FP-Zernike/, as well as a webserver http://43.138.38.101.
## Packages && deploy
    In order to run FP-Zernike, three main packages are required, that is, 'pymol', 'numpy' and 'numba'. For convenience, an environment file('FP_Zernike.yaml') is provided so that users can deploy FP-Zernike directly using conda, the command line is as follows:
        conda env create -f FP_Zernike.yaml
        pip install scikit-learn
    In addtion, users need to create a folder by command line as follows:
        mkdir tempoutput
    And, open permissions for some tools:
        chmod 777 tool/gmconvert
        chmod 777 tool/extractFeaturePoints30
        chmod 777 tool/extractFeaturePoints50

## Usage
### 1. Calculate Zernike descriptors
    Given the query structure[path: example/4mu3A.pdb], the user can execute the following command to get its ATOM-Zernike:
        python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source atom -mode atom -o test_4mu3A_ATOM.pkl
    Here, 'test_4m2qA_ATOM.pkl' include the information of Zernike descriptor of 'example/4m2qA.pdb'. 
    The user can execute the following command to get its PM-Zernike:
        python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source pymol -mode mesh -o test_4mu3A_PM.pkl
    The user can execute the following command to get its PS-Zernike:
        python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source pymol -mode surface -o test_4mu3A_PS.pkl
    The user can execute the following command to get its GM-Zernike:
        python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source gmconvert -mode mesh -o test_4mu3A_GM.pkl
    In addition, users can enter the following command to obtain detailed instructions.
        python mainZernikeCalculation.py -h
    The above commands in 'test_mainZernikeCalculation.sh', user can run 'bash test_mainZernikeCalculation.sh' to test this function.
### 2. Build a database of descriptors
    Given the customized dataset[only include *.pdb file, path: example], the user can execute the following command to get its ATOM-Zernike database:
        python mainMakeZernikeDataBase.py -i example -t 4 -mode ATOM
    Here, '4' is the thread number. A new folder[path: exampleAtomZernikeDescriptor] containing all descriptors is automatically generated. Similarly, descriptors for other patterns can also be built, and these command lines are in 'test_mainMakeZernikeDataBase.sh', user can run 'bash test_mainMakeZernikeDataBase.sh' to test this function.
    In addition, users can enter the following command to obtain detailed instructions.
        python mainMakeZernikeDataBase.py -h

### 3. Ultra-fast structure retrieval
    Given the query structure, the user can perform ultra-fast structure retrieval. For example, given a protein chain[Path:example/4mu1A.pdb], users can preform sturcture retrieval by command as follows.
        python mainRetrieve.py -q example/4mu1A.pdb -dm PM -otn 10 -op testPMRtriRes.txt
    Here, 'PM' is the mode of retrieve, '10' is the number of output structures, 'testPMRtriRes.txt' stores all the output.
    In addition, users can enter the following command to obtain detailed instructions.
        python mainRetrieve.py -h
    There are more commands in 'test_mainRetrieve.sh', and users can test this function by executing 'bash test_mainRetrieve.sh'.
    We constructed a demo containing 590685 structures, at this scale, our retrieval system only takes 4 ∼ 9 seconds to complete a retrieval. To perform the retrieval on the demo dataset, users need to through the following link to download about demo database (Name:PMAllSingleChainZernike.pkl), and need to put     'PMAllSingleChainZernike.pkl' copy to the folder 'database'.
        link: https://pan.baidu.com/s/1NnQsFiVyatQ2p-k7oj-giA 
        Extraction code: 208z

### 4. Measure the similarity between the two structures
    When the descriptors of the structure are calculated by FP-Zernike, the user can initially measure the similarity based on the Euclidian distance between the descriptors. Specifically, given the descriptors of two structures, A[Path: example/4mr5A.pkl] and B[Path: example/4mr6A.pkl], the user can perform the following command to obtain the Euclidean distance between them.
        python mainZernikeDist.py example/4mr5A.pkl example/4mr6A.pkl
### 5. Measure the similarity based on ReOmokage
    Here, we improve the Omokage and propose 'ReOmokage' to measure the similarity between structures. Given structures A[Path: example/4mu3A.pdb] and B[Path: example/4mu1A.pdb], the user can enter the following command to calculate the similarity between the structures:
        python mainReOmokage.py -pdb1 example/4mu3A.pdb -pdb2 example/4mu1A.pdb

### cite
Qi, J., Feng, C., Shi, Y., Yang, J., Zhang, F., Li, G., & Han, R. (2024). FP-Zernike: an open-source structural database construction toolkit for fast structure retrieval. Genomics, Proteomics & Bioinformatics, qzae007.

    

        





