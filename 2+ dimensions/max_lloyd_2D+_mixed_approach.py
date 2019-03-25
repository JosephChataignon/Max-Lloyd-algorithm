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
    for k in xrange(nbdimensions):
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

## function studied
def f(position):
    return gaussian(position)
    #return uniform(position)

## distribution studied
def random_distrib():
    #return [random.uniform(-1,1),random.uniform(-1,1)]
    return [random.gauss(0,1),random.gauss(0,1)]

## squared euclidian distance between x and y for any number of dimensions
def sqdistance(x,y):
    d = 0.0
    for i in xrange(len(x)):
        d += (x[i]-y[i])**2
    return d

## returns the index of the region of which point is part of
def isinregion(point,regions):
    for k in xrange(len(regions)): # for each region k
        isin = True
        for l in xrange(len(regions)): #for each region l distinct from k
            if k != l:
                if np.dot(regions[k,l],np.append(point,1.0)) > 0: #if point out of bounds
                    isin = False
        if isin:
            return k
    return -1

## returns the mean squared error on the region between -5 and +5 for each axis
# As for the centroid function that comes after, it is an approximation only
# calculated on points of a grid in a limited area
def MSE(regions, germs, griddimensions):
    error = 0.
    total_proba = 0. # cumulated probabilities in region k, initialized as 0
    # Here is the same loop structure as in the centroid function
    indexes = [0 for x in xrange(len(griddimensions))]  # initialize the loops indexes
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
        index_coord = [float(indexes[x])*10./griddimensions[x]-5. for x in xrange(len(indexes))] # correspondance between indexes and coordinate system
        regioncenter = germs[isinregion(index_coord,regions)]
        proba_density = f(index_coord)
        total_proba += proba_density
        error += proba_density * sqdistance( regioncenter , index_coord )
    if total_proba != 0: # should be equal to 1
        error /= total_proba
    else:
        print "ERROR: total_proba = 0"
        coord = np.zeros(len(griddimensions))
    return error

## This function computes the centroid of the k-th region,
# based on an estimation in each point of a grid. It is thus an approximation.
# In particular it only takes into account a 10-width zone centered on zero 
# and not the parts of the region that have coordinates <-5 or >5
def centroid(regions,k,griddimensions):
    coord = np.zeros(len(griddimensions))
    total_proba = 0. # cumulated probabilities in region k, initialized as 0
    # Here comes a particular loop structure, because the number of dimensions
    # (ie the number of imbricated loops we must use) is variable
    indexes = [0 for x in xrange(len(griddimensions))]  # initialize the loops indexes
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
        index_coord = [float(indexes[x])*10./griddimensions[x]-5. for x in xrange(len(indexes))] # correspondance between indexes and coordinate system
        if isinregion(index_coord,regions) == k :
            proba_density = f(index_coord)
            total_proba += proba_density
            for x in xrange(len(griddimensions)):
                coord[x] += proba_density * index_coord[x]
    if total_proba != 0:
        coord /= total_proba
    else:
        print "total_proba = 0"
        coord = np.zeros(len(griddimensions))
    return coord

## This function adjusts the delimitations of the regions based on their center
# it builds a Voronoi diagram from an array containing a germ for each region
def adjust_regions(germs):
    nbdimensions = len(germs[0])    # number of dimensions
    nbregions = len(germs)          # number of regions (or number of germs)
    regions = np.array( [ [ [0.0]*(nbdimensions+1) ]*nbregions ]*nbregions )
    for i in xrange(nbregions):
        for j in xrange(nbregions):
            if i == j:
                pass
            else:
                b = 0.0
                for k in xrange(nbdimensions):
                    regions[i,j,k] = 2 * ( germs[j,k] - germs[i,k] )
                    b += germs[i,k]**2 - germs[j,k]**2
                regions[i,j,-1] = b
                if b != 0:
                    regions[i,j,:] = regions[i,j,:] / np.abs(b)
    return regions

## This function is the actual implementation of Max-Lloyd algorithm
# regions delimiters are a set of variables [a1,a2,...an,b] , which are 
# parameters of the equation a1*x1 + a2*x2 + ... + an*xn + b <= 0
# griddimensions indicates the number of points on each axis that are used for
# estimating the centroid of a region
def maxlloyd(germs,griddimensions,iterations):
    c = 0                           # counts the number of iterations
    nbdimensions = len(germs[0])    # number of dimensions
    nbregions = len(germs)          # number of regions
    regions = np.array( [ [ [0.0]*(nbdimensions+1) ]*nbregions ]*nbregions )
    while c < iterations:
        c = c+1
        if c%2 == 1:
            # adjust regions
            regions = adjust_regions(germs)
        else:
            # adjust germs
            for i in xrange(len(germs)):
                germs[i] = centroid(regions,i,griddimensions)
    return germs,regions


## function for visualising the regions, 2D only, for a number of regions <= 8
def displayregions(griddimensions,regions):
    plt.figure(2)
    colors = ['b', 'c', 'y', 'm', 'r', 'g', 'w', 'k']
    for i in xrange(griddimensions[0]):
        for j in xrange(griddimensions[1]):
            i2 = float(i)*10.0/griddimensions[0]-5 ; j2 = float(j)*10.0/griddimensions[1]-5; point = np.array([i2,j2])
            plt.plot(i2,j2,marker='o',color=colors[isinregion(point,regions)])
    plt.show()

## Test for displayregions function, only with regions from initial germs
def test_displayregions():
    germs = np.array([[1.0,1.0],[-1.0,2.0],[0.0,0.0],[2.0,1.0],[1.0,-1.0],[-1.0,-2.0],[3.0,-2.0],[3.0,1.0]])
    regions = adjust_regions(germs)
    print regions,"\n"
    displayregions([100,100],regions)
#test_displayregions() # uncomment to test


## Test for maxlloyd function
def test_maxlloyd():
    germs = np.array([[1.,1.],[-1.,2.] ,[0.,0.],[2.,1.],[1.,-1.],[-1.,-2.],[3.,-2.],[3.,1.]])
    grid_dim = [60,60]
    iterations = 20
    germs,regions = maxlloyd(germs,grid_dim,iterations)
    displayregions(grid_dim,regions)
    MSE(regions,germs,grid_dim)
    return germs,regions
#test_maxlloyd() # uncomment to test

## This function measures the MSE depending on the number of regions
def measure_error(iterations,griddimensions,number_of_regions,number_of_repeats):
    error = 0.
    for i in xrange(number_of_repeats):
        germs = []
        for k in xrange(number_of_regions):
            germs.append([random.uniform(-3,3) for z in xrange(len(griddimensions))])
        germs = np.array(germs)
        germs,regions = maxlloyd(germs,griddimensions,iterations)
        error += MSE(regions,germs,griddimensions)
    error /= number_of_repeats
    return error
#print measure_error(10,[15,15],4,3)

## This function measures the evolution of MSE over increasing number of regions
def error_evolution(iterations,griddimensions,number_of_repeats,encoding_size):
    error_values = []
    for n in xrange(1,encoding_size):
        error_values.append(measure_error(iterations,griddimensions,2**n,number_of_repeats))
    return error_values
#print error_evolution(10,[15,15,15],3)

## This function displays the evolution of MSE for several dimensions
# gridresolution is the number of points on one axis of the grid
# iterations is the number of iterations on maxlloyd function
# number of repeats is the number of times calculations are repeated to measure average error
# n_dimensions is the number of dimensions for which error evolution is computed
# encoding size is the maximal number of bits for the quantizer output
# ie for n<encoding_size, the number of regions is 2**n
def display_error(iterations,gridresolution,number_of_repeats,n_dimensions,encoding_size):
    plt.figure(1)
    for d in xrange(1,n_dimensions):
        print "dimension",d
        griddimensions = [gridresolution for i in xrange(d)]
        error = error_evolution(iterations,griddimensions,number_of_repeats,encoding_size)
        plt.plot(error,label='%d dimensions'%d)
    plt.legend(loc='upper right')
    plt.show()
display_error(iterations=20,gridresolution=100,number_of_repeats=10,n_dimensions=5,encoding_size=8)











