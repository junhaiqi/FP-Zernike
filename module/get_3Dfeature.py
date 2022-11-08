# coding=utf-8
import os
from Bio.Data.SCOPData import protein_letters_3to1
import numpy as np
# import get_position_ofInterface as gpo
# from double_check_extint import getInterfacePosition as gIp


def get_atom_coordinate(pdb_file):
    '''
    description: given a pdb file, output heavy atom coordinates about residues.
    '''

    all_coords = {}
    # loop pdb lines, construct all atom seqs
    with open(pdb_file, 'r') as fr:
        last_chain, last_aa_num, aa_seq_coord, chain_seq_coord = '', '', [], []
        for line in fr:
            if line.startswith('ATOM'):
                # parse ATOM line
                aa = line[17:20]  # example: 'ALA'
                if aa in protein_letters_3to1.keys():
                    chain = line[21]  # example: 'A'
                    aa_num = line[22:27].strip()  # example: '1'
                    # find heavy atom coordinate
                    # line[13] indicates atom type, example: 'C'.
                    atom_coord = []
                    if line[13] != 'H':
                        temp_atom_coord = line[27:54].strip().split(' ')
                        try:
                            atom_coord = tuple([float(item) for item in temp_atom_coord if
                                                item != ''])
                        # throw away bad coordinate, e.g., '49.680-100.424  42.761'
                        except:
                            print('throw away 3D coordinate %s' %
                                  str(line[27:54].strip()))
                    if atom_coord != []:
                        # find new chain
                        if chain != last_chain:
                            chain_seq_coord.append(aa_seq_coord)
                            all_coords[last_chain] = chain_seq_coord
                            chain_seq_coord = []
                            aa_seq_coord = []
                            # aa_seq_coord.append(atom_coord)
                            last_chain = chain
                        # find new amino acid
                        if aa_num == last_aa_num:
                            aa_seq_coord.append(atom_coord)
                            # last_aa_num, last_chain = aa_num, chain
                        else:
                            chain_seq_coord.append(aa_seq_coord)
                            aa_seq_coord = [atom_coord]
                            last_aa_num, last_chain = aa_num, chain
                            # last_chain = chain
                else:
                    print(
                        'Warning: pdb file include non-protein structure and this structure will be omitted!')
        # loop end
        chain_seq_coord.append(aa_seq_coord)
        all_coords[last_chain] = chain_seq_coord
        # remove redundant elements in the key value of dictionary.
        for key in all_coords:
            all_coords[key].pop(0)
        # remove empty key value of dictionary.
        all_coords.pop('')
        return all_coords


def EuclideanDistance(tuple: tuple, tuple_list: list):
    '''
    description: a function for calculating the distance list between a tuple and a tuple list.

    args:
        tuple: a tuple, example: (1,2,3)
        tuple_list: a tuple list, example: [(1,2,3),(1,2,4),(1,2,1)]

    output: a distance list. for example, tuple: (1,2,3), tuple list: [(1,2,3),(1,2,4),(1,2,1)]
            output: [0,1,2]
    '''

    norm = np.linalg.norm
    array = np.array
    dist_list = norm(array(tuple) - array(tuple_list), axis=1)

    return min(dist_list)


def getFeature_Chains(pdb_file, chain1, chain2, chain3, pair_cutoff=10):
    '''
    description: given a pdb file, output interface residues position of all chains.

    args:
        pdb_file: input file
        chain1: 'A' chain, antigen chain.
        chain2: 'H' chain, heavy chain.
        chain3: 'L' chain, light chain.
        pair_cutoff: a distance cutoff for defining point pair.

    output:
        distChain_dict  ---distance dictionary, it indicates Residue-to-chain distance.
        For example, {'A': array([33.71482738,...,), 'H': ...}, here, '33.71482738' 
        indicates the distance from the first residue in 'A' to 'H' and 'L' is '33.71482738'.

        point_pair_list  ---point pair list. For example, '(12, 53, 9.691331435876084)'
        indicates the distance from the 12th residue in 'A' to 53th residue in 'H' and 'L' is
        9.691331435876084[less than defalut cutoff=10].
    '''

    # get atom coordinate of residues in chains.
    atom_coord_dict = get_atom_coordinate(pdb_file)
    # get the coordinate list about 'A' chain.
    chain_coord1 = atom_coord_dict[chain1]
    # get the coordinate list about 'HL' chain, that is consider HL two chains as one.
    chain_coord2 = atom_coord_dict[chain2] + atom_coord_dict[chain3]
    # a length value for spliting 'H' and 'L' information.
    split_HL_len = len(atom_coord_dict[chain2])
    # a distance matric that indicates the distance from the residue in chain 'A' to the residue in chain 'H' or 'L'.
    dist_matric = np.zeros((len(chain_coord1), len(chain_coord2)))
    # interface contact pair list.
    # point_pair_list1 ---load contact pair between 'A' and 'H'. point_pair_list2 ---load contact pair between 'A' and 'L'.
    point_pair_list1, point_pair_list2 = [], []
    # Global variable to temporary variable.
    _append = point_pair_list1.append
    _append1 = point_pair_list2.append
    _min = min
    _min1 = np.min

    # for loop for getting distance matric and contact pair between 'A' and 'HL'.
    for i in range(0, len(chain_coord1)):
        for j in range(0, len(chain_coord2)):

            dist_list = [EuclideanDistance(
                item, chain_coord2[j]) for item in chain_coord1[i]]
            min_dist = _min(dist_list)
            dist_matric[i][j] = min_dist

            if min_dist < pair_cutoff:
                if j < split_HL_len:
                    _append((i, j, min_dist))
                else:
                    _append1((i, j-split_HL_len, min_dist))

    # creat a distance dictionary for three chains.
    distChain_dict = {}
    # get the distance vector about chain 'A'. This vector indicates the distance from residue in 'A' to 'HL'.
    distChain_dict[chain1] = _min1(dist_matric, axis=1)
    # split 'H' and 'L' by the length of chain 'H'[The first 'split_HL_len' positions belong to 'H'].
    distChain_dict[chain2] = _min1(dist_matric, axis=0)[0:split_HL_len]
    # [The rest positions belong to 'L']
    distChain_dict[chain3] = _min1(dist_matric, axis=0)[split_HL_len:]
    # creat a triplet dictionary for loading contact pair. example: (1,2,2.3), the third position in '(1,2,2.3)' indicate distance.
    contact_pair_dict = {}
    # name
    key_1 = chain1 + '_' + chain2
    key_2 = chain1 + '_' + chain3
    # contact pair between chain 'A' and chain 'H'.
    contact_pair_dict[key_1] = point_pair_list1
    # contact pair between chain 'A' and chain 'L'.
    contact_pair_dict[key_2] = point_pair_list2

    # get interface dictionary.
    interface_dict = {}
    interface_dict[chain1] = [1 if item <
                              pair_cutoff else 0 for item in distChain_dict[chain1]]
    interface_dict[chain2] = [1 if item <
                              pair_cutoff else 0 for item in distChain_dict[chain2]]
    interface_dict[chain3] = [1 if item <
                              pair_cutoff else 0 for item in distChain_dict[chain3]]

    # return, get distance dictionary, contact pair dictionary and interface dictionary.
    return distChain_dict, contact_pair_dict, interface_dict


def get_splitcount(k, list):

    if k == 0:
        return (0, list[0])
    else:
        SplitCount = 0
        for i in range(0, k):
            SplitCount += list[i]

        return (SplitCount, SplitCount+list[k])


def get_range_all_k(list):
    range_list = [get_splitcount(k, list) for k in range(0, len(list))]
    return range_list


def split_antigen_dist_list(dist_list, antigen_len_list):

    dist_list_split = []
    for i in range(0, len(antigen_len_list)):
        split_index = get_splitcount(i, antigen_len_list)
        # print(split_index)
        dist_list_split.append(dist_list[split_index[0]:split_index[1]])

    return dist_list_split


def upgrade_contact_info(contact_list, k, antigen_len_list):
    '''
    description:
                contact_list: a contact list
                k: indicate the antigen chain
                antigen_len_list: a list for loading the length of each antigen chain
    '''
    # print('**',len())
    pre_count = get_splitcount(k, antigen_len_list)[0]
    new_contact_list = [(item[0]-pre_count, item[1], item[2])
                        for item in contact_list]

    return new_contact_list


def GetFeatureChainsListVersion(pdb_file, antigen_list, chain2, chain3, pair_cutoff=10):
    '''
    description: given a pdb file, output interface residues position of all chains.

    args:
        pdb_file: input file
        antigen_list: antigen chains[antigen has multiple chains.].
        chain2: 'H' chain, heavy chain.
        chain3: 'L' chain, light chain.
        pair_cutoff: a distance cutoff for defining point pair.

    output:
        distChain_dict  ---distance dictionary, it indicates Residue-to-chain distance.
        For example, {'A': array([33.71482738,...,), 'H': ...}, here, '33.71482738' 
        indicates the distance from the first residue in 'A' to 'H' and 'L' is '33.71482738'.

        point_pair_list  ---point pair list. For example, '(12, 53, 9.691331435876084)'
        indicates the distance from the 12th residue in 'A' to 53th residue in 'H' and 'L' is
        9.691331435876084[less than defalut cutoff=10].
    '''

    # a list for loading the length of each antigen chain.
    antigen_chain_len_list = []
    # get atom coordinate of residues in chains.
    atom_coord_dict = get_atom_coordinate(pdb_file)
    # get the coordinate list about antigen chains.
    chain_coord1 = []
    for item in antigen_list:
        antigen_chain_len_list.append(len(atom_coord_dict[item]))
        chain_coord1 += atom_coord_dict[item]
    print('chain_A_len:', len(chain_coord1))
    # get the coordinate list about 'HL' chain, that is consider HL two chains as one.
    chain_coord2 = atom_coord_dict[chain2] + atom_coord_dict[chain3]
    # a length value for spliting 'H' and 'L' information.
    split_HL_len = len(atom_coord_dict[chain2])
    # a distance matric that indicates the distance from the residue in chain 'A' to the residue in chain 'H' or 'L'.
    dist_matric = np.zeros((len(chain_coord1), len(chain_coord2)))
    # interface contact pair list.
    # point_pair_list1 ---load contact pair between 'A' and 'H'. point_pair_list2 ---load contact pair between 'A' and 'L'.
    point_pair_list1, point_pair_list2 = [], []
    # Global variable to temporary variable.
    _append = point_pair_list1.append
    _append1 = point_pair_list2.append
    _min = min
    _min1 = np.min

    # for loop for getting distance matric and contact pair between 'A' and 'HL'.
    for i in range(0, len(chain_coord1)):
        for j in range(0, len(chain_coord2)):

            dist_list = [EuclideanDistance(
                item, chain_coord2[j]) for item in chain_coord1[i]]
            min_dist = _min(dist_list)
            dist_matric[i][j] = min_dist

            if min_dist < pair_cutoff:
                if j < split_HL_len:
                    _append((i, j, min_dist))
                else:
                    _append1((i, j-split_HL_len, min_dist))

    # creat a distance dictionary for three chains.
    distChain_dict = {}
    # get the distance vector about chain 'A'. This vector indicates the distance from residue in 'A' to 'HL'.
    distChain_dict_antigen = _min1(dist_matric, axis=1)
    antigen_dist_list = split_antigen_dist_list(
        distChain_dict_antigen, antigen_chain_len_list)
    for k in range(0, len(antigen_list)):
        distChain_dict[antigen_list[k]] = antigen_dist_list[k]
    # split 'H' and 'L' by the length of chain 'H'[The first 'split_HL_len' positions belong to 'H'].
    distChain_dict[chain2] = _min1(dist_matric, axis=0)[0:split_HL_len]
    # [The rest positions belong to 'L']
    distChain_dict[chain3] = _min1(dist_matric, axis=0)[split_HL_len:]
    # creat a triplet dictionary for loading contact pair. example: (1,2,2.3), the third position in '(1,2,2.3)' indicate distance.
    contact_pair_dict = {}
    range_list = get_range_all_k(antigen_chain_len_list)
    # get the key names of 'contact_pair_dict'
    for _range_index in range(0, len(range_list)):
        pair_name1 = antigen_list[_range_index] + '_H'
        pair_name2 = antigen_list[_range_index] + '_L'
        contact_pair_dict[pair_name1], contact_pair_dict[pair_name2] = [], []

    # upgrade contact pair information of 'A_H'
    for pair in point_pair_list1:
        for _range_index in range(0, len(range_list)):
            if pair[0] >= range_list[_range_index][0] and pair[0] < range_list[_range_index][1]:
                pair_name = antigen_list[_range_index] + '_H'
                contact_pair_dict[pair_name].append(
                    (pair[0]-range_list[_range_index][0], pair[1], pair[2]))

    # upgrade contact pair information of 'A_L'
    for pair in point_pair_list2:
        for _range_index in range(0, len(range_list)):
            if pair[0] >= range_list[_range_index][0] and pair[0] < range_list[_range_index][0]:
                antigen_list[_range_index] + '_L'
                contact_pair_dict[pair_name].append(
                    (pair[0]-range_list[_range_index][0], pair[1], pair[2]))

    # get interface dictionary.
    interface_dict = {}
    for k in range(0, len(antigen_list)):
        interface_dict[antigen_list[k]] = [1 if item <
                                           pair_cutoff else 0 for item in antigen_dist_list[k]]
    # interface_dict[chain1] = [1 if item <
    #                           pair_cutoff else 0 for item in distChain_dict[chain1]]
    interface_dict[chain2] = [1 if item <
                              pair_cutoff else 0 for item in distChain_dict[chain2]]
    interface_dict[chain3] = [1 if item <
                              pair_cutoff else 0 for item in distChain_dict[chain3]]

    # return, get distance dictionary, contact pair dictionary and interface dictionary.
    return distChain_dict, contact_pair_dict, interface_dict


if __name__ == "__main__":

    get_atom_coordinate('../ranked_0.pdb')
    # current_file_path = os.path.dirname(os.path.abspath('.'))
    # print(current_file_path)
