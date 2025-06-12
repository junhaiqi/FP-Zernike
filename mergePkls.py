import pickle
import os
import sys

if not os.path.exists('database'): 
    os.makedirs(f'database')   
    
def extractInfoFromPklfile(pkl_file):
    
    f = open(pkl_file, 'rb')
    data = pickle.load(f)
    f.close()
    print( data )     
          
def merge_pkl_to_dict(database_dir_path, output_file):
    
    pkl_files = [ database_dir_path + '/' + item for item in os.listdir(database_dir_path) ]
    
    merged_dict = {}
    
    for file_path in pkl_files:
        if not os.path.exists(file_path):
            continue

        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        
        if isinstance(data, dict):
            merged_dict.update(data)
        else:
            key = os.path.splitext(os.path.basename(file_path))[0]  # 去掉扩展名
            merged_dict[key] = data
    
    with open(output_file, 'wb') as f:
        pickle.dump(merged_dict, f)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} $Your_Zernike_Descriptor_Dir_Path $Output_Pkl_File")
        exit(-1)
    else:
        merge_pkl_to_dict(database_dir_path = sys.argv[1], output_file = sys.argv[2]')
    
