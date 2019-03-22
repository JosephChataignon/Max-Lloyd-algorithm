# -*- coding: utf-8 -*-
#
# Max-Lloyd algorithm for finding an optimal quantizer
# in dimension 2 or more
#
# This is a mixed approach to try to get results more efficiently than with the
# numerical approach, and extended to higher dimensions


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

## squared euclidian distance between x and y
def sqdistance(x,y):
    d = 0.0
    for i in xrange(len(x)):
        d += (x[i]-y[i])**2
    return d

## returns the index of the region of which point is part of
def isinregion(point,regions):
    for k in xrange(len(regions)):
        isin = True
        for l in xrange(len(regions)):
            if k != l:
                if np.dot(regions[k,l],np.append(point,1.0)) > 0: #if point out of bounds
                    isin = False
        if isin:
            return k
    return -1

## This function computes the centroid of the k-th region,
# based on an estimation in each point of a grid
def centroid(regions,k,griddimensions):
    coord = np.zeros(len(griddimensions))
    total_proba = 0.0 # cumulated probabilities in the k-th region
    
    ## Here comes a particular loop structure, because we don't know in advance 
    # how many dimensions there are (ie how many imbricated loops we must use)
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
        ## Actual loop content
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

## This function adjusts the delimitations of each region based on their center
# it builds a Voronoi diagram with the array of germs
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
        print "iteration:",c
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


## Test for maxlloyd function
def test_maxlloyd():
    germs = np.array([[1.,1.],[-1.,2.] ,[0.,0.],[2.,1.],[1.,-1.],[-1.,-2.],[3.,-2.],[3.,1.]])
    grid_dim = [30,30]
    iterations = 10
    germs,regions = maxlloyd(germs,grid_dim,iterations)
    displayregions(grid_dim,regions)
    return germs,regions

test_maxlloyd()




def test_displayregions():
    #germs = np.array([[1.0,1.0],[-1.0,2.0] ,[0.0,0.0],[2.0,1.0]])
    germs = np.array([[1.0,1.0],[-1.0,2.0],[0.0,0.0],[2.0,1.0],[1.0,-1.0],[-1.0,-2.0],[3.0,-2.0],[3.0,1.0]])
    regions = adjust_regions(germs)
    print regions,"\n"
    displayregions([30,30],regions)
    
#test_displayregions()






















# the structure for variable number of loops has been inspired by this example:
# http://python.jpvweb.com/python/mesrecettespython/doku.php?id=boucles_imbriquees