# FP-Zernike: An open-source structural database construction toolkit for fast and accurate structure retrieval

The main functions of FP-Zernike as follows:

1. Given a protein(or RNA) structure(format: *.pdb), FP-Zernike can calculate its Zernike descriptors[4 types, namely PM-Zernike, PS-Zernike, ATOM-Zernike and GM-Zernike]. Based on the results of the experiment, we recommend the user to calculate PM-Zernike, PS-Zernike or ATOM-Zernike.

2. Given a customized dataset[only include *.pdb file], FP-Zernike can calculate all Zernike descriptors about it.

3. Given a query structure[single-chain], structure retrieval is done on a dataset containing 590685 protein chains[run time: ~4s-9s].

4. Based on the Euclidean distance of descriptors, the similarity of structures is preliminarily evaluated. When the Euclidean distance between descriptors is less than 2, it means that there is strong similarity.

5. Based on 'Romokage', the similarity between structures is measured.

FP-Zernike is an open-source toolkit, whose source code and related data are accessible at https://github.com/junhaiqi/FP-Zernike/.

## Packages && deploy

Installation Option 1: In order to run FP-Zernike, three main packages are required, that is, 'pymol', 'numpy' and 'numba'. For convenience, an environment file('FP_Zernike.yaml') is provided so that users can deploy FP-Zernike directly using conda, the command line is as follows:
    
```bash
    conda env create -f PGAR_Zernike.yaml
    conda activate PGAR_Zernike
    pip install scikit-learn
```
        
Installation Option 2: The command line is as follows:
    
```bash
    conda create -n FP-Zernike python==3.7
    conda activate FP-Zernike
    conda install numpy==1.20 numba==0.55.2
    conda install conda-forge::pymol-open-source
    pip install scipy tqdm
```
        
In addtion, users need to create a folder by command line as follows:

```bash
    mkdir tempoutput
```
    
And, open permissions for some tools:
    
```bash
    chmod 777 tool/gmconvert
    chmod 777 tool/extractFeaturePoints30
    chmod 777 tool/extractFeaturePoints50
 ```

## Usage
### 1. Calculate Zernike descriptors

Given the query structure[path: example/4mu3A.pdb], the user can execute the following command to get its ATOM-Zernike:

 ```bash
    python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source atom -mode atom -o test_4mu3A_ATOM.pkl
```

Here, 'test_4m2qA_ATOM.pkl' include the information of Zernike descriptor of 'example/4m2qA.pdb'. 

The user can execute the following command to get its PM-Zernike:

 ```bash
    python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source pymol -mode mesh -o test_4mu3A_PM.pkl
```

The user can execute the following command to get its PS-Zernike:

 ```bash
    python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source pymol -mode surface -o test_4mu3A_PS.pkl
 ```

The user can execute the following command to get its GM-Zernike:

 ```bash
    python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source gmconvert -mode mesh -o test_4mu3A_GM.pkl
 ```

In addition, users can enter the following command to obtain detailed instructions.

 ```bash
    python mainZernikeCalculation.py -h
 ```

The above commands in 'test_mainZernikeCalculation.sh', user can run 'bash test_mainZernikeCalculation.sh' to test this function.

### 2. Build a database of descriptors

Given the customized dataset[only include *.pdb file, path: example], the user can execute the following command to get its ATOM-Zernike database:

 ```bash
    python mainMakeZernikeDataBase.py -i example -t 4 -mode ATOM
```

Here, '4' is the thread number. A new folder[path: exampleAtomZernikeDescriptor] containing all descriptors is automatically generated. Similarly, descriptors for other patterns can also be built, and these command lines are in 'test_mainMakeZernikeDataBase.sh', user can run 'bash test_mainMakeZernikeDataBase.sh' to test this function.

In addition, users can enter the following command to obtain detailed instructions.

 ```bash
    python mainMakeZernikeDataBase.py -h
```

### 3. Ultra-fast structure retrieval

Given the query structure, the user can perform ultra-fast structure retrieval. For example, given a protein chain[Path:example/4mu1A.pdb], users can preform sturcture retrieval by command as follows.

 ```bash
    python mainRetrieve.py -q example/4mu1A.pdb -dm PM -otn 10 -op testPMRtriRes.txt
```

Here, 'PM' is the mode of retrieve, '10' is the number of output structures, 'testPMRtriRes.txt' stores all the output.

In addition, users can enter the following command to obtain detailed instructions.

 ```bash
    python mainRetrieve.py -h
```

There are more commands in 'test_mainRetrieve.sh', and users can test this function by executing 'bash test_mainRetrieve.sh'.

We constructed a demo containing 590685 structures, at this scale, our retrieval system only takes 4 ∼ 9 seconds to complete a retrieval. To perform the retrieval on the demo dataset, users need to through the following link to download about demo database (Name:PMAllSingleChainZernike.pkl), and need to copy 'PMAllSingleChainZernike.pkl' to the folder 'database'.

 ```bash
    Qi, Junhai (2025). FP-Zernike Descriptor Database. figshare. Dataset. https://doi.org/10.6084/m9.figshare.29304539.v1
```

If you have generated a descriptor database (assuming the path is `example` and the mode is `PM`) and want to perform a retrieve based on it, you can use the following command line configuration:

```
    python mergePkls.py example database/PMAllSingleChainZernike.pkl
```

Here, the naming rule is database/$modeAllSingleChainZernike.pkl.

### 4. Measure the similarity between the two structures
When the descriptors of the structure are calculated by FP-Zernike, the user can initially measure the similarity based on the Euclidian distance between the descriptors. Specifically, given the descriptors of two structures, A[Path: example/4mr5A.pkl] and B[Path: example/4mr6A.pkl], the user can perform the following command to obtain the Euclidean distance between them.

 ```bash
    python mainZernikeDist.py example/4mr5A.pkl example/4mr6A.pkl
```

### 5. Measure the similarity based on ReOmokage

Here, we improve the Omokage and propose 'ReOmokage' to measure the similarity between structures. Given structures A[Path: example/4mu3A.pdb] and B[Path: example/4mu1A.pdb], the user can enter the following command to calculate the similarity between the structures:

 ```bash
    python mainReOmokage.py -pdb1 example/4mu3A.pdb -pdb2 example/4mu1A.pdb
```

## Cite
Junhai Qi, Chenjie Feng, Yulin Shi, Jianyi Yang, Fa Zhang, Guojun Li, Renmin Han, FP-Zernike: An Open-source Structural Database Construction Toolkit for Fast Structure Retrieval, Genomics, Proteomics & Bioinformatics, Volume 22, Issue 1, February 2024, qzae007, https://doi.org/10.1093/gpbjnl/qzae007

    

        





