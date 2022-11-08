python mainZernikeCalculation.py -i example/4m2qA.pdb -file_source atom -mode atom -o test_4m2qA_ATOM.pkl
python mainZernikeCalculation.py -i example/4m2qA.pdb -file_source pymol -mode mesh -o test_4m2qA_PM.pkl
python mainZernikeCalculation.py -i example/4m2qA.pdb -file_source pymol -mode surface -o test_4m2qA_PS.pkl
python mainZernikeCalculation.py -i example/4m2qA.pdb -file_source gmconvert -mode GM -o test_4m2qA_GM.pkl
