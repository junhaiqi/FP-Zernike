
'''
@Author  :   qijunhai
@Email   :   qijunhai@mail.sdu.edu.cn
@description: This python script for calculating 3D zernike descriptor by PGAR-Zernike algorithm. 
              The difinition and calculation metails of 3D zernike descriptor can be found 
              in following article.
              
              "3D zernike descriptors for content based shape retrieval,2003"             
'''
import os
import numpy as np
import math
from scipy.special import (factorial, comb as nchoosek)
from tqdm import tqdm
import numba as nb
nb.config.NUMBA_DEFAULT_NUM_THREADS = 8
import sys
import argparse

if 'module' not in sys.path:
    sys.path.append('module')

import pickle
# import use_tool
import get_protein_representation as gpr
import time


def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True,
                        help='input file, it is a pdb file.')
    parser.add_argument('-o', type=str, required=True,
                        help='output file that used to store zernike descriptor(format: pkl).')
    parser.add_argument('-file_source', type=str, required=True,
                        help='indicate the source of protein representation(gmconvert, atom or pymol).')
    parser.add_argument('-mode', type=str, required=True,
                        help='protein representation type. The following options are available: surface, mesh, atom, sphere(only -file_source pymol)')
    parser.add_argument('-gmmCount', type=int, required=False, default=50, 
                        help='a parameter to gmconvert that controls the resolution of the resulting surface.')

    args = parser.parse_args()
    return args

@nb.njit(parallel=True, nogil=True)
def getMoments(index_list, diff=0.01):

    # index_list = np.array(index_list)
    # time_start = time.time()
    
    moments = np.zeros((21, 21, 21))
    
    print('##############################################################')
    print('&&&&&&&&&&&&&&Start Calculating Geometric Moment&&&&&&&&&&&&&&')
    for r in nb.prange(0, 21):
        print('Now r =', r, ', when r = 20, the geometric moment is calculated.')
        for s in nb.prange(0, 21):
            for t in nb.prange(0, 21):
                if r + s + t < 21:
                    
                    divide = 1 / (r+1) * 1/(s+1) * 1/(t+1)
                    moment = 0
                    for item in index_list:
                        x_step = ((item[0] + 1)*diff - 1)**(r + 1) - (item[0]*diff - 1)**(r+1)
                        y_step = ((item[1] + 1)*diff - 1)**(s + 1) - (item[1]*diff - 1)**(s+1)
                        z_step = ((item[2] + 1)*diff - 1)**(t + 1) - (item[2]*diff - 1)**(t+1)
                        
                        moment += x_step*y_step*z_step*divide
                        
                    moments[r][s][t] = moment
                    
    # time_end = time.time()

    print('##############################################################')
    print('&&&&&&&&&&&&&&&End Calculating Geometric Moment&&&&&&&&&&&&&&&')
    
    return moments


class cal_moments_rst():

    def __init__(self, surface_file_name, mode='wrl', file_source='pymol'):
        '''
        Initial parameters function, only need one parameter("surface_file_name").
        '''

        self.surface_file_name = surface_file_name  # a '*.wrl' file
        self.mode = mode  # This parameter for adding some function in future.
        # This parameter for indicating the source of protein representation file.
        self.file_source = file_source

    def get_atom_coordinate(self, mode = 'noheavy'):
        
        all_atom_coord = []
        with open(self.surface_file_name, 'r') as fr:
            for line in fr:
                
                if line.startswith('TER'):
                    break
                
                elif line.startswith('ATOM'):
                    if mode == 'heavy':
                        if line[77] != 'H':
                            temp_atom_coord = line[27:54].strip().split(' ')
                            try:
                                atom_coord = tuple([float(item) for item in temp_atom_coord if
                                                    item != ''])

                                all_atom_coord.append(atom_coord)
                            # throw away bad coordinate, e.g., '49.680-100.424  42.761'
                            except:
                                print('throw away 3D coordinate %s' %
                                        str(line[27:54].strip()))
                    else:
                        temp_atom_coord = line[27:54].strip().split(' ')
                        try:
                            atom_coord = tuple([float(item) for item in temp_atom_coord if
                                                item != ''])

                            all_atom_coord.append(atom_coord)
                        # throw away bad coordinate, e.g., '49.680-100.424  42.761'
                        except:
                            print('throw away 3D coordinate %s' %
                                    str(line[27:54].strip()))
                        
        if all_atom_coord == []:
            print('Warning: point set is empty!')
        return all_atom_coord

    # def getPointCloudCoordinate(self):
    #     pointList = []
    #     # surface_file_name is the point cloud file.
    #     with open(self.surface_file_name, 'r') as fr:
    #         pass
        
    #     return pointList

    # @nb.jit()
    def get_points(self):
        # start = time.time()
        '''
        description: This simple function for getting 3D point cloud from a '*.wrl' file.
                     '*.wrl' file is obtained by 'gmconvert' or 'pymol'(two tools).
        '''
        if self.file_source == 'atom':
            points_list = self.get_atom_coordinate()
            return np.array(points_list)

        # elif self.file_source == 'pointCloud':
        #     points_list = self.getPointCloudCoordinate()
        #     return np.array(points_list)

        elif self.file_source == 'hybrid':
            points_list1 = self.get_atom_coordinate()

            pre_name = self.surface_file_name.split('/')[-1]
            pre_name = pre_name.split('.')[0]
            
            gmfileList = os.listdir(path='tempoutput')
            if '%s_ng100_surface.wrl'%pre_name not in gmfileList:
                use_tool.get_pdbSurface(self.surface_file_name)
            
            with open('tempoutput/%s_ng100_surface.wrl'%pre_name, 'r') as surface_file:
                lines = surface_file.readlines()
                index_1 = 0
                index_2 = 0
                
                str_index_1 = ' point[\n'
                str_index_2 = ']}\n'
                
                for i in range(0, len(lines)):
                    # print(lines[i])
                    if lines[i] == str_index_1:
                        # print(lines[i])
                        index_1 = i + 1

                    if lines[i] == str_index_2:
                        index_2 = i
                        break
                    
                points_str = lines[index_1:index_2]

                # extrat 3D points from lines.
                points_list2 = []
                for item in points_str:
                    item = item.split('\n')
                    item = item[0].split(',')
                    item = item[0].split(' ')
                    temp_list = [float(point_str) for point_str in item]
                    points_list2.append(temp_list)
                    
            return np.array(points_list1), np.array(points_list2)
        
        else:

            with open(self.surface_file_name, 'r') as surface_file:
                lines = surface_file.readlines()
                index_1 = 0
                index_2 = 0
                
                if self.file_source == 'gmconvert':
                    str_index_1 = ' point[\n'
                    str_index_2 = ']}\n'
                elif self.file_source == 'pymol':
                    str_index_1 = '   point [\n'
                    str_index_2 = '   ]\n'

                for i in range(0, len(lines)):
                    # print(lines[i])
                    if lines[i] == str_index_1:
                        # print(lines[i])
                        index_1 = i + 1

                    if lines[i] == str_index_2:
                        index_2 = i
                        break

                # From '*.wrl' file gets lines with 3D points.
                points_str = lines[index_1:index_2]

                # extrat 3D points from lines.
                points_list = []
                
                _append = points_list.append
                for item in points_str:
                    item = item.split('\n')
                    item = item[0].split(',')
                    item = item[0].split(' ')
                    temp_list = [float(point_str) for point_str in item]
                    _append(temp_list)

                # return a array with 3D points
                return np.array(points_list)

    def scale_to_unit_ball(self, mode='wrl'):
        '''
        description: First, get 3D points from a input file(*.wrl). Then, the 3D points
                     are scaled into a unit ball(x^2+y^2+z^2<=1).
        '''

        if self.file_source != 'hybrid':
            
        # get 3D points from a input file(*.wrl).

            print('Extracting information for 3D points...')
            start_getPoints = time.time()
            points = self.get_points()
            end_Points = time.time()
            print('The time taken to extract the information of three-dimensional points is: %ss.' \
                  %(end_Points - start_getPoints))
            
            if len(points) == 0:
                print('File format error(Atom points or surface points)! Calculation failed!')
                exit()
            # weighted mean to get a point as centroid.
            centroid = sum(points) / len(points)
            # The centroid is translated to origin, and other points translate the corresponding position.
            points = points - centroid
            _norm = np.linalg.norm
            scale_index =  max([item for item in _norm(points, ord = 2, axis = 1, keepdims = True)])
        
            scale_points = points / (scale_index + 1)
            
            return scale_points
        
        else:
            points1 = self.get_points()[0]
            points2 = self.get_points()[1]
            
            centroid1 = sum(points1)/len(points1)
            points1 = points1 - centroid1
            
            centroid2 = sum(points2)/len(points2)
            points2 = points2 - centroid2
            
            scale_index1 = np.linalg.norm(points1[0])
            for item in points1:
                if np.linalg.norm(item) > scale_index1:
                    scale_index1 = np.linalg.norm(item)

            # The 3D points are scaled into a unit ball.
            scale_points1 = points1 / (scale_index1 + 1)
            
            scale_index2 = np.linalg.norm(points2[0])
            for item in points2:
                if np.linalg.norm(item) > scale_index2:
                    scale_index2 = np.linalg.norm(item)

            # The 3D points are scaled into a unit ball.
            scale_points2 = points2 / (scale_index2 + 1)
            
            scale_points = np.array(list(scale_points1) + list(scale_points2))
            
            # return
            return scale_points

    def get_grid_index(self, coordinate, N=200):
        '''
        description：There is a cube, its center is origin and consists of 
                     200 * 200 * 200 small cubes, the length of small cube 
                     is 2/200, so we can know the length of big cube is 2.
                     This function do this: given a 3D coordinate, determine 
                     which small cube this coordinate is in. Note: this input 
                     coordinate in unit ball(center is origin). The range of 
                     cube is [-1,1] for x, y and z coordinate.          
        '''

        diff = 2 / N  # the length of small cube.
        x = 0
        y = 0
        z = 0
        
        x = int(coordinate[0] / diff) + int(N/2)
        y = int(coordinate[1] / diff) + int(N/2)
        z = int(coordinate[2] / diff) + int(N/2)
        
        # return
        return [x, y, z]

    def get_function_about_cube(self, N=200):
        '''
        description: First, get the 3D points are scaled. Then, determine which small cube 
                     these 3D points are in. 
        '''

        time_start = time.time()
        # print('Number of points: ', len(points))
        print('Start placing the point in the unit sphere...')
        points = self.scale_to_unit_ball()  # get 3D points
        print('Number of points:', len(points), '.')
        # print('points count: ', len(points))
        grid_valueIsone_list = []  # a list consist of the small cubes that include 3D point
        time_end = time.time()
        print('The time it takes to place the point in the unit sphere is: %ss.'%(time_end - time_start))
        
        # _append = grid_valueIsone_list.append
        
        print('Start building voxels...')
    
        start_valueOne = time.time()
        
        grid_valueIsone_list = [self.get_grid_index(item) for item in points]
            
        end_valueOne = time.time()
        
        print('The time of the voxelization is: %ss.'%(end_valueOne - start_valueOne))
            
        # return
        return grid_valueIsone_list

    def getIndexList(self, order=20, N=200, mode='surface'):
        '''
        description: calculate geometric moment(a type of integral). Please refer 
                     to paper given in this script description for more metails.

        input: 'order' indicate the order of 3D zernike moment.
                N=200 is a default value.
                mode='surface' is a default.
        output: a list consists of geometric moment for 'rst'.
        '''

        # print('start CUBE...\n')
        start_getCube = time.time()
        index_list = self.get_function_about_cube()
        end_getCube = time.time()
        # print('CUBE Time: %f'%(end_getCube - start_getCube))
        
        return np.array(index_list)


class Zernike():
    '''
    description: given moments that calculated by "cal_moments_rst()", to calculate
                 3D zernike descriptors.
    '''

    def __init__(self, moments):

        self.binomial_dict = {}
        self.factorial_dict = {}
        self.q_klv_dict = {}
        self.clm_dict = {}
        self.moments = moments

    def get_nlms(self, order):
        '''
        description: a function for getting a list that consists of triple that
                     certain conditions are satisfied. 
        '''

        nlms = []
        for n in range(order+1):
            for l in range(n+1):
                for m in range(l+1):
                    if (n-l) % 2 == 0:
                        nlms.append([n, l, m])

        return np.array(nlms)

    def find_same_nl(self, nlm_list):
        '''
        description: Group a list that consists of nlm triples. The element of 
                     a group have same 'n' and 'l'.

        example: input [[1,2,3],[1,2,4],[3,4,1],[2,3,5]], output
                       [[0,1],2,3]        
        '''

        # a list for getting nlm triples index that have been grouped.
        index_list = []
        output_list = []  # output
        for i in range(0, len(nlm_list)-1):
            if i not in index_list:
                temp_list = []
                temp_list.append(i)
                index_list.append(i)
                for j in range(i+1, len(nlm_list)):
                    if self.find_same_count(nlm_list[i], nlm_list[j]) > 1:
                        temp_list.append(j)
                        index_list.append(j)
                output_list.append(temp_list)

        return output_list

    def find_same_count(self, list_1, list_2):

        # a simple function for helping "find_same_nl(self,nlm_list)"
        t = 0
        if list_1[0] == list_2[0] and list_1[1] == list_2[1]:
            t = 2

        return t

    def binomial(self, n, k):
        '''
        description: a function for calculating combination(C(n,k))
        '''

        if (n, k) in self.binomial_dict.keys():
            output = self.binomial_dict[(n, k)]
        else:
            output = nchoosek(n, k)
            self.binomial_dict[(n, k)] = output

        return output

    def fact(self, n):
        '''
        description: a function for calculating factorial(n!)
        '''
        if n in self.factorial_dict.keys():
            output = self.factorial_dict[n]
        else:
            output = factorial(n)
            self.factorial_dict[n] = output

        return output

    '''
    The following functions are for calculation “chi_nlm_rst”, The calculation process designs a 
    large number of cycles. Please refer to peper:
             "3D zernike descriptors for content based shape retrieval,2003"
    '''

    def chi_nlm_rst(self, n, l, m, mode):

        output = self.c_l_m(l, m) * (2.**(-m)) + 0j
        k = (n-l)//2
        output = output * self.sum1(n, l, m, k, mode)
        if mode == 'omega':
            output = np.conj(output) * 0.75 / np.pi

        return output

    def c_l_m(self, l, m):

        if (l, m) in self.clm_dict.keys():
            output = self.clm_dict[(l, m)]
        else:
            output = (2*l + 1) * self.fact(l+m) * self.fact(l-m)
            output = math.sqrt(output)
            output /= self.fact(l)

            self.clm_dict[(l, m)] = output

        return output

    def sum1(self, n, l, m, k, mode):

        output = 0 + 0j
        for nu in range(k+1):
            output += self.Q_klv(k, l, nu) * self.sum2(n, l, m, k, nu, mode)

        return output

    def Q_klv(self, k, l, nu):

        if (k, l, nu) in self.q_klv_dict.keys():
            output = self.q_klv_dict[(k, l, nu)]

        else:
            output = (-1)**(k+nu)
            output /= 2.**(2*k)
            output *= math.sqrt((2*l + 4*k + 3)/3)
            output *= self.binomial(2*k, k)
            output *= self.binomial(k, nu)
            output *= self.binomial(2*(k+l+nu)+1, 2*k)
            output /= self.binomial(k+l+nu, k)
            self.q_klv_dict[(k, l, nu)] = output

        return output

    def sum2(self, n, l, m, k, nu, mode):

        output = 0 + 0j
        for alpha in range(nu+1):
            output += self.binomial(nu, alpha) * \
                self.sum3(n, l, m, k, nu, alpha, mode)

        return output

    def sum3(self, n, l, m, k, nu, alpha, mode):

        output = 0 + 0j
        for beta in range(nu-alpha+1):
            output += self.binomial(nu-alpha, beta) * \
                self.sum4(n, l, m, k, nu, alpha, beta, mode)

        return output

    def sum4(self, n, l, m, k, nu, alpha, beta, mode):

        output = 0 + 0j
        for u in range(m+1):
            output += ((-1.)**(m-u)) * self.binomial(m, u) * ((0 + 1j)**u) * self.sum5(
                n, l, m, nu, alpha, beta, u, mode)

        output = output.real-(output.imag)*(0+1j)

        return output

    def sum5(self, n, l, m, nu, alpha, beta, u, mode):

        output = 0 + 0j
        for mu in range(((l-m)//2)+1):
            temp = ((-1.)**mu) * (2.**(-2*mu))
            temp *= self.binomial(l, mu)
            temp *= self.binomial(l-mu, m+mu)
            temp *= self.sum6(n, l, m, nu, alpha, beta, u, mu, mode)
            output += temp

        return output

    def sum6(self, n, l, m, nu, alpha, beta, u, mu, mode):

        output = 0 + 0j
        for v in range(mu+1):
            r = 2*(v+alpha)+u
            s = 2*(mu-v+beta)+m-u
            t = 2*(nu-alpha-beta-mu)+l-m
            temp = self.binomial(mu, v)
            if mode == 'omega':
                temp *= self.moments[r, s, t]
                output += temp

        return output

    def calculate_chi(self, order, mode='omega'):

        self.q_klv_dict = {}
        self.clm_dict = {}
        nlms = self.get_nlms(order)
        chi_d = {}

        for triplet in nlms:
            n, l, m = triplet
            chi = self.chi_nlm_rst(n, l, m, mode)
            chi_d[(n, l, m)] = chi

        return chi_d

def inv_zscore(inv_arrary):
    
    '''
    description: a function for data standardization.
    '''

    mean = sum(inv_arrary) / len(inv_arrary)
    omega = (sum((inv_arrary - mean) ** 2) / len(inv_arrary)) ** 0.5
    inv_zscore_array = (inv_arrary - mean) / omega
    return inv_zscore_array

# defalut order = 20

def get_zernike(file_source, input, gmmCount, mode, output):
    
    '''
    description: Connect the above functions in series to realize 
                 the functions: input pdb file and output zernike
                 descriptors.
    '''

    zernike_start = time.time()
    order = 20
    # get structure surface by "gmconvert"
    if file_source == 'gmconvert':
        pre_name = input.split('/')[-1]
        pre_name = pre_name.split('.')[0]
        # gmfileList = os.listdir(path='tempoutput')
        # if '%s_ng50_surface.wrl'%pre_name not in gmfileList:
        start_mesh = time.time()
        print('Generating gmconvert protein representation file...')
        Gmmcount = gmmCount
        use_tool.get_pdbSurface(input, GMMCount=Gmmcount)
        end_mesh = time.time()
        
        print('The time taken to generate the protein representation file is: %ss.'%(end_mesh - start_mesh))
        
        # calculate geometric moment.
        momentsClass = cal_moments_rst(surface_file_name='tempoutput/%s_ng%d_surface.wrl' % (pre_name, Gmmcount),
                               file_source='gmconvert')
        
        indexList = momentsClass.getIndexList()
        
        start_time = time.time()
        moments = getMoments(indexList)
        end_time = time.time()
        print('The calculation time of geometric moment is: %ss.'%(end_time - start_time))

    elif file_source == 'pymol':
        pre_name = input.split('/')[-1]
        pre_name = pre_name.split('.')[0]
        start_mesh = time.time()
        print('Generating pymol protein representation file...')
        gpr.main(input, 'tempoutput/%s_%s.wrl' %
                 (pre_name, mode), mode=mode)
        end_mesh = time.time()
        
        print('The time taken to generate the protein representation file is: %ss.'%(end_mesh - start_mesh))
        
        momentsClass = cal_moments_rst(surface_file_name='tempoutput/%s_%s.wrl' % (pre_name, mode),
                               mode=mode, file_source=file_source)
        
        indexList = momentsClass.getIndexList()
        
        start_time = time.time()
        moments = getMoments(indexList)
        end_time = time.time()
        print('The calculation time of geometric moment is: %ss.'%(end_time - start_time))

    elif file_source == 'atom':
        
        pre_name = input.split('/')[-1]
        pre_name = pre_name.split('.')[0]
        
        momentsClass = cal_moments_rst(surface_file_name=input,
                               mode=mode, file_source=file_source)
        
        indexList = momentsClass.getIndexList()
        
        start_time = time.time()
        moments = getMoments(indexList)
        end_time = time.time()
        print('The calculation time of geometric moment is: %ss.'%(end_time - start_time))
    
    elif file_source == 'hybrid':
        pre_name = input.split('/')[-1]
        pre_name = pre_name.split('.')[0]
        # gmfileList = os.listdir(path='tempoutput')
       
        momentsClass = cal_moments_rst(surface_file_name = input,
                               file_source='hybrid')
        
        indexList = momentsClass.getIndexList()
        
        start_time = time.time()
        moments = getMoments(indexList)
        end_time = time.time()
        print('The calculation time of geometric moment is: %ss.'%(end_time - start_time))

    # calculate chi dictionary
    ZernikeMoments = Zernike(moments)
    s = ZernikeMoments.get_nlms(order)
    same = ZernikeMoments.find_same_nl(s)
    discriptor_dict = ZernikeMoments.calculate_chi(order)

    nl_list = []
    for item in same:
        temp_list = []
        for item_1 in item:
            temp_list.append(discriptor_dict[tuple(s[item_1])])
        nl_list.append(temp_list)

    # calculate 3D zernike descriptors
    inv_list = []
    for item in nl_list:
        inv = 0
        for item_1 in item:
            inv += abs(item_1) ** 2
        inv_list.append(math.sqrt(inv))

    # z-score 3D zernike descriptors
    inv_list = inv_zscore(np.array(inv_list))
    inv_dict = {input.split('/')[-1].split('.')[0]: inv_list}
    file = open(output, 'ab')
    pickle.dump(inv_dict, file, -1)
    file.close()
    
    zernike_end = time.time()
    print('Zernike feature vector has been written into %s. Total time: %ss.' % (output, zernike_end - zernike_start))

    return inv_list

def main():
    args = initialization_parameters()
    get_zernike(file_source = args.file_source, input = args.i, gmmCount = args.gmmCount, mode = args.mode, output = args.o)



if __name__ == '__main__':

    # get_zernike()
    main()
    