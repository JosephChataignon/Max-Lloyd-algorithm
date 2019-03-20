# -*- coding: utf-8 -*-
#
# Max-Lloyd algorithm for finding the optimal quantizer
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
    return np.exp(-0.5*position**2) / ((2.0*np.pi)**(len(position)*0.5))

# function studied
def f(position):
    return gaussian(position)
    #return uniform(position)

# distribution studied
def random_distrib():
    #return [random.uniform(-1,1),random.uniform(-1,1)]
    return [random.gauss(0,1),random.gauss(0,1)]

# squared euclidian distance between x and y
def sqdistance(x,y):
    d = 0.0
    for i in xrange(len(x)):
        d = d + (x[i]-y[i])**2
    return d

# k is the index of the germ to adjust
def centroid(grid,k):
    Cx = 0.0; Cy = 0.0; total_proba = 0.0
    for i in xrange(len(grid)):
        for j in xrange(len(grid[i])):
            if grid[i,j] == k:
                proba_density = f( float(i)*10.0/float(len(grid))-5.0 , float(j)*10.0/float(len(grid[0]))-5.0 )
                total_proba = total_proba + proba_density
                Cx = Cx + proba_density * (float(i)*10.0/float(len(grid))-5.0)
                Cy = Cy + proba_density * (float(j)*10.0/float(len(grid[0]))-5.0)
    if total_proba != 0:
        Cx = Cx / total_proba
        Cy = Cy / total_proba
    else:
        Cx = 0; Cy = 0
    return [Cx,Cy]

# germs is an array containing all the germs
def adjust_regions(germs):
    nbdimensions = len(germs[0]) # number of dimensions
    nbregions = len(germs) # number of regions (or number of germs)
    regions = np.array( [ [ [0.0]*(nbdimensions+1) ]*nbregions ]*nbregions )
    for i in xrange(nbregions):
        for j in xrange(nbregions):
            if i == j:
                #to complete
                print "not implemented yet"
            else:
                b = 0.0
                for k in xrange(nbdimensions):
                    regions[i,j,k] = 2 * ( germs[j,k] - germs[i,k] )
                    b = b + germs[i,k]**2 - germs[j,k]**2
                regions[i,j,-1] = b
                if b != 0:
                    regions[i,j,:] = regions[i,j,:] / np.abs(b)
    return regions


# regions delimiters are a set of variables [a1,a2,...an,b] , which are 
# parameters of the equation a1*x1 + a2*x2 + ... + an*xn + b <= 0
def maxlloyd(germs,griddimensions,iterations):
    c = 0                           # counts the number of iterations
    nbdimensions = len(germs[0])    # number of dimensions
    nbregions = len(germs)          # number of regions
    regions = np.array( [ [ [0.0]*(nbdimensions+1) ]*nbregions ]*nbregions )
    while c < iterations:
        c = c+1
        print c
        if c%2 == 1:
            # adjust regions
            grid = adjust_regions(germs)
        else:
            # adjust germs
            for i in xrange(len(germs)):
                germs[i] = centroid(regions,i)
    return germs,grid


# Test of maxlloyd function
def test_maxlloyd():
    #germs = [[1.0,1.0],[-1.0,2.0],[0.0,0.0],[2.0,1.0]]
    germs = np.array([[1.0,1.0],[-1.0,2.0],[0.0,0.0],[2.0,1.0],[1.0,-1.0],[-1.0,-2.0],[3.0,-2.0],[3.0,1.0]])
    grid_dim = [100,100]
    iterations = 100
    germs,grid = maxlloyd(germs,grid_dim,iterations)
    plt.figure(2)
    colors = ['b', 'c', 'y', 'm', 'r', 'g', 'k', 'w']
    for i in xrange(len(grid)):
        for j in xrange(len(grid[0])):
            plt.plot(float(i)*10.0/15.0-5.0,float(j)*10.0/15.0-5.0,marker='o',color=colors[grid[i,j]])
    plt.show()
    return germs,grid

#test_maxlloyd()


# returns the index of the region of which point is part of
def isinregion(point,regions):
    for k in xrange(len(regions)):
        isin = True
        for l in xrange(len(regions)):
            if k != l:
                if np.dot(regions[k,l],np.append(point,1.0)) > 0: #if point out of bounds
                    isin = False
        if isin:
            return k
    return 4
    

# test function for visualising the regions of the grid after adjust_grid
def displayregions():
    germs = np.array([[1.0,1.0],[-1.0,2.0] ,[0.0,0.0],[2.0,1.0]])
    regions = adjust_regions(germs)
    print regions,"\n"
    plt.figure(2)
    for i in xrange(-5,5):
        for j in xrange(-5,5):
            i2 = float(i) ; j2 = float(j); point = np.array([i2,j2])
            if isinregion(point,regions) == 0:
                plt.plot(i2,j2,'ro')
            elif isinregion(point,regions) == 1:
                plt.plot(i2,j2,'go')
            elif isinregion(point,regions) == 2:
                plt.plot(i2,j2,'yo')
            elif isinregion(point,regions) == 3:
                plt.plot(i2,j2,'bo')
            else:
                plt.plot(i2,j2,'ko')
    plt.show()

displayregions()

























