# -*- coding: utf-8 -*-
#
# Max-Lloyd algorithm for finding an optimal quantizer
# in dimension 2 or more
#
# This is a mixed approach to try to get results more efficiently than with the
# numerical approach, and extended to higher dimensions
#
# the structure for variable number of loops has been inspired by this example:
# http://python.jpvweb.com/python/mesrecettespython/doku.php?id=boucles_imbriquees


import math
import random
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



def uniform(position):
    nbdimensions = len(position)
    in_distribution = True
    for k in range(nbdimensions):
        if position[k]<-1 or position[k]>1:
            in_distribution = False
    if in_distribution:
        return 2.0**(-nbdimensions)
    else:
        return 0.0

def gaussian(position):
    p = 1.
    for coord in position:
        p *= np.exp(-0.5*coord**2) / ((2.0*np.pi)**(len(position)*0.5))
    return p

def f(position):
    '''function studied'''
    return gaussian(position)
    #return uniform(position)

def random_distrib():
    '''distribution studied'''
    #return [random.uniform(-1,1),random.uniform(-1,1)]
    return [random.gauss(0,1),random.gauss(0,1)]

def sqdistance(x,y):
    '''squared euclidian distance between x and y for any number of dimensions'''
    d = 0.0
    for i in range(len(x)):
        d += (x[i]-y[i])**2
    return d

def isinregion(point,regions):
    '''returns the index of the region of which point is part of'''
    for k in range(len(regions)): # for each region k
        isin = True
        for l in range(len(regions)): #for each region l distinct from k
            if k != l:
                if np.dot(regions[k,l],np.append(point,1.0)) > 0: #if point out of bounds
                    isin = False
        if isin:
            return k
    return -1

def MSE(regions, germs, griddimensions):
    '''
        Returns the mean squared error on the region between -5 and +5 for each 
        axis.
        As for the centroid function that comes after, it is an approximation 
        only calculated on points of a grid in a limited area.
    '''
    error = 0.
    total_proba = 0. # cumulated probabilities in region k, initialized as 0
    # Here is the same loop structure as in the centroid function, see below
    indexes = [0 for x in range(len(griddimensions))]  # initialize the loops indexes
    indexes[-1] -= 1                                    # before 1st incrementation
    indmax = len(griddimensions)-1  # index of the last loop (most internal one)
    endloop = False                 # boolean for ending condition of the loop
    while True:
        for x in range(indmax, -1, -1): # reverse loop of loops number (from internal to external)
            indexes[x] += 1
            if indexes[x] >= griddimensions[x]: # if loop x's index is out of bounds
                if x == 0: # if the first loop's end has been reached
                    endloop = True
                    break
                indexes[x] = 0
            else: break # index[x] is inside bounds, break and execute loop content
        if endloop: break
        # Actual loop content
        # index_coord is a point of the space - correspondance between indexes and coordinate system:
        index_coord = [float(indexes[x])*10./griddimensions[x]-5. for x in range(len(indexes))]
        # find the centroid of the region this point is belonging to
        regioncenter = germs[isinregion(index_coord,regions)]
        # get the probability density on that point
        proba_density = f(index_coord)
        # add it to the total of probabilities
        total_proba += proba_density
        # augment error with the appropriate value
        error += proba_density * sqdistance( regioncenter , index_coord )
    if total_proba != 0: # total_proba should be equal to 1
        error /= total_proba
    else:
        print("ERROR: total_proba = 0")
        coord = np.zeros(len(griddimensions))
    return error

def centroid(regions,k,griddimensions):
    '''
        This function computes the centroid of the k-th region, based on an 
        estimation in each point of a grid. It is thus an approximation. In 
        particular it only takes into account a 10-width zone centered on zero 
        and not the parts of the region that have coordinates <-5 or >5
    '''
    coord = np.zeros(len(griddimensions))
    total_proba = 0. # cumulated probabilities in region k, initialized as 0
    # Here comes a particular loop structure, because the number of dimensions
    # (ie the number of imbricated loops we must use) is variable
    indexes = [0 for x in range(len(griddimensions))]  # initialize the loops indexes
    indexes[-1] -= 1                                    # before 1st incrementation
    indmax = len(griddimensions)-1  # index of the last loop (most internal one)
    endloop = False                 # boolean for ending condition of the loop
    while True:
        for x in range(indmax, -1, -1): # reverse loop of loops number (from internal to external)
            indexes[x] += 1
            if indexes[x] >= griddimensions[x]: # if loop x's index is out of bounds
                if x == 0: # if the first loop's end has been reached
                    endloop = True
                    break
                indexes[x] = 0
            else: break # index[x] is inside bounds, break and execute loop content
        if endloop: break
        # Actual loop content
        index_coord = [float(indexes[x])*10./griddimensions[x]-5. for x in range(len(indexes))] # correspondance between indexes and coordinate system
        if isinregion(index_coord,regions) == k :
            proba_density = f(index_coord)
            total_proba += proba_density
            for x in range(len(griddimensions)):
                coord[x] += proba_density * index_coord[x]
    if total_proba != 0:
        coord /= total_proba
    else:
        print("total_proba = 0")
        coord = np.zeros(len(griddimensions))
    return coord

def adjust_regions(germs):
    '''
        This function adjusts the delimitations of the regions based on their
        center. It builds a Voronoi diagram from an array containing a germ for 
        each region.
    '''
    nbdimensions = len(germs[0])    # number of dimensions
    nbregions = len(germs)          # number of regions (or number of germs)
    regions = np.array( [ [ [0.0]*(nbdimensions+1) ]*nbregions ]*nbregions )
    for i in range(nbregions):
        for j in range(nbregions):
            if i == j:
                pass
            else:
                b = 0.0
                for k in range(nbdimensions):
                    regions[i,j,k] = 2 * ( germs[j,k] - germs[i,k] )
                    b += germs[i,k]**2 - germs[j,k]**2
                regions[i,j,-1] = b
                if b != 0:
                    regions[i,j,:] = regions[i,j,:] / np.abs(b)
    return regions

def maxlloyd(germs,griddimensions,iterations):
    '''
        Actual implementation of Max-Lloyd algorithm.
        Regions delimiters are a set of variables [a1,a2,...an,b] , which are 
        parameters of the equation a1*x1 + a2*x2 + ... + an*xn + b <= 0
        griddimensions indicates the number of points on each axis that are used for
        estimating the centroid of a region.
    '''
    c = 0                           # counts the number of iterations
    nbdimensions = len(germs[0])    # number of dimensions
    nbregions = len(germs)          # number of regions
    regions = np.array( [ [ [0.0]*(nbdimensions+1) ]*nbregions ]*nbregions )
    print("number of regions:",len(germs))
    while c < iterations:
        c = c+1
        if c%2 == 1:
            # adjust regions
            regions = adjust_regions(germs)
        else:
            # adjust germs
            for i in range(len(germs)):
                germs[i] = centroid(regions,i,griddimensions)
    return germs,regions

def performance():
    '''generate figures with maxlloyd function for the paper'''
    #germs = np.array([[0.,-1.],[0.,1.]])#2
    #germs = np.array([[0.,-1.],[-0.75,0.25],[0.75,0.25]])#3
    #germs = np.array([[0.5,-0.5],[-0.75,-0.75],[-0.25,0.]])#3bis
    #germs = np.array([[0.,1.],[-1.,0.],[1.,0.],[0.,-1.]])#4
    #germs = np.array([[-0.5,-0.5],[-0.5,0.5],[0.5,0.5],[0.5,-0.5]])#4bis
    #germs = np.array([[-0.75,0.],[-0.25,0.8],[0.4,0.4],[0.7,-0.2]])#4ter
    #germs = np.array([[-0.75,0.6],[-0.2,-0.3],[0.25,-0.63],[0.87,0.12]])#4quad
    #germs = np.array([[-0.75,0.6],[-0.2,-0.3],[0.25,-0.63],[0.87,0.12],[-0.36,0.46],[0.83,-0.35],[-0.52,0.84] ])#7
    
    regionsnumber = 16
    germs = np.reshape([np.random.random(regionsnumber*2)*2-1],(regionsnumber,2))
    
    germs, regions = maxlloyd(germs,[100,100],30)
    m = MSE(regions,germs,[1000,1000])
    print(m)
    displayregions([30,30],regions)

def displayregions(griddimensions,regions):
    '''visualizes the regions, in 2D only, for a number of regions <= 8'''
    plt.figure(2)
    colors = ['b', 'c', 'y', 'm', 'r', 'g', 'w', 'k']
    for i in range(griddimensions[0]):
        for j in range(griddimensions[1]):
            i2 = float(i)*10.0/griddimensions[0]-5 ; j2 = float(j)*10.0/griddimensions[1]-5; point = np.array([i2,j2])
            plt.plot(i2,j2,marker='o',color=colors[isinregion(point,regions)])
    plt.show()

performance()

def test_displayregions():
    '''Test for displayregions function, only with regions from initial germs'''
    germs = np.array([[1.0,1.0],[-1.0,2.0],[0.0,0.0],[2.0,1.0],[1.0,-1.0],[-1.0,-2.0],[3.0,-2.0],[3.0,1.0]])
    regions = adjust_regions(germs)
    print(regions,"\n")
    displayregions([100,100],regions)
#test_displayregions() # uncomment to test


def test_maxlloyd():
    '''Test for maxlloyd function'''
    germs = np.array([[1.,1.],[-1.,2.] ,[0.,0.],[2.,1.],[1.,-1.],[-1.,-2.],[3.,-2.],[3.,1.]])
    grid_dim = [60,60]
    iterations = 20
    germs,regions = maxlloyd(germs,grid_dim,iterations)
    displayregions(grid_dim,regions)
    MSE(regions,germs,grid_dim)
    return germs,regions
#test_maxlloyd() # uncomment to test

def measure_error(iterations,griddimensions,number_of_regions,number_of_repeats):
    '''measures the MSE depending on the number of regions'''
    error = 0.
    for i in range(number_of_repeats):
        germs = []
        for k in range(number_of_regions):
            germs.append([random.uniform(-3,3) for z in range(len(griddimensions))])
        germs = np.array(germs)
        germs,regions = maxlloyd(germs,griddimensions,iterations)
        error += MSE(regions,germs,griddimensions)
    error /= number_of_repeats
    return error
#print measure_error(10,[15,15],4,3)

def error_evolution(iterations,griddimensions,number_of_repeats,encoding_size):
    '''measures the evolution of MSE over increasing number of regions'''
    error_values = []
    for n in range(1,encoding_size+1):
        error_values.append(measure_error(iterations,griddimensions,n**len(griddimensions),number_of_repeats))
    return error_values
#print error_evolution(10,[15,15,15],3)

def display_error(iterations,gridresolution,number_of_repeats,n_dimensions,encoding_size):
    '''
        This function displays the evolution of MSE for several dimensions.
        gridresolution is the number of points on one axis of the grid.
        iterations is the number of iterations on maxlloyd function.
        number of repeats is the number of times calculations are repeated to 
        measure average error.
        n_dimensions is the number of dimensions for which error evolution is 
        computed.
        encoding_size is the maximal number of regions per axis,
        ie for n <= encoding_size, the number of regions is n**dimensions
    '''
    plt.figure(1)
    for d in range(1,n_dimensions):
        print("dimension",d)
        griddimensions = [gridresolution for i in range(d)]
        if d == 3: # last 3 are not computed due to very high computation time needed
            error = error_evolution(iterations,griddimensions,number_of_repeats,encoding_size-3)
        else:
            error = error_evolution(iterations,griddimensions,number_of_repeats,encoding_size)
        plt.plot(error,label='%d dimensions'%d)
    plt.legend(loc='upper right')
    plt.show()
#display_error(iterations=10,gridresolution=50,number_of_repeats=1,n_dimensions=4,encoding_size=6)











