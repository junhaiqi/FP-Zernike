# Test the build descriptor database. Here, 'example' is a folder path that include structures(format: pdb), '1' is the number of thread.
python mainMakeZernikeDataBase.py -i example -t 1
#Here, mode is the type of Zernike descriptor.
python mainMakeZernikeDataBase.py -i example -t 1 -mode ATOM
python mainMakeZernikeDataBase.py -i example -t 1 -mode PM
python mainMakeZernikeDataBase.py -i example -t 1 -mode PS
python mainMakeZernikeDataBase.py -i example -t 1 -mode GM
