python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source atom -mode atom -o test_4mu3A_ATOM.pkl
python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source pymol -mode mesh -o test_4mu3A_PM.pkl
python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source pymol -mode surface -o test_4mu3A_PS.pkl
python mainZernikeCalculation.py -i example/4mu3A.pdb -file_source gmconvert -mode surface -o test_4mu3A_GM.pkl
