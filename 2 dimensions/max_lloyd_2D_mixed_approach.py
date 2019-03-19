
# Max-Lloyd algorithm for finding the optimal quantizer
# in dimension 2

import math
import random
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate



def uniform(x,y):
    if x<=1 and x>=-1 and y<=1 and y>=-1:
        return 0.25
    else:
        return 0.0

def gaussian(x,y):
    return math.exp((-(float(x)**2.0)-(float(y)**2.0))/2.0)

# function studied
def f(x,y):
    #return gaussian(x,y)
    return uniform(x,y)

# distribution studied
def random_distrib():
    #return [random.uniform(-1,1),random.uniform(-1,1)]
    return [random.gauss(0,1),random.gauss(0,1)]

# squared distance between 2 points
def sqdistance(point1,point2):
    return (float(point1[0])-float(point2[0]))**2 + (float(point1[1])-float(point2[1]))**2


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

#germs is an array containing all the germs
def adjust_grid(germs,grid_dim):
    grid = np.array([[0]*grid_dim[0]]*grid_dim[1])
    for i in xrange(len(grid)):
        for j in xrange(len(grid[0])):
            point = [float(i)*10.0/float(len(grid))-5.0 , float(j)*10.0/float(len(grid[0]))-5.0]
            nearest = 0
            for k in xrange(1,len(germs)):
                if sqdistance(point,germs[k]) < sqdistance(point,germs[nearest]):
                    nearest = k
            grid[i,j] = nearest
    return grid


# grid is a regular grid which covers the plane around the point (0,0). The
# integer contained at each point is the index of the nearest germ from this point
# germs is an array containing the germ for each region, on the model (x,y)
# error_threshold is the threshold to reach for the algorithm to stop
def maxlloyd(germs,grid_dim):
    print "test fonction  maxlloyd"
    c = 0 # counts the number of executions of the loop
    grid = np.array([[0]*grid_dim[0]]*grid_dim[1])
    while c < 100:
        c = c+1
        print c
        if c%2 == 1:
            # adjust grid
            grid = adjust_grid(germs,grid_dim)
        else:
            # adjust germs
            for i in xrange(len(germs)):
                germs[i] = centroid(grid,i)
    return germs,grid


# Test of maxlloyd function
def test_maxlloyd():
    #germs = [[1.0,1.0],[-1.0,2.0],[0.0,0.0],[2.0,1.0]]
    germs = [[1.0,1.0],[-1.0,2.0],[0.0,0.0],[2.0,1.0],[1.0,-1.0],[-1.0,-2.0],[3.0,-2.0],[3.0,1.0]]
    grid_dim = [100,100]
    germs,grid = maxlloyd(germs,grid_dim)
    plt.figure(2)
    colors = ['b', 'c', 'y', 'm', 'r', 'g', 'k', 'w']
    for i in xrange(len(grid)):
        for j in xrange(len(grid[0])):
            plt.plot(float(i)*10.0/15.0-5.0,float(j)*10.0/15.0-5.0,marker='o',color=colors[grid[i,j]])
    plt.show()
    return germs,grid

test_maxlloyd()



# test function for visualising the regions of the grid after adjust_grid
def displayregions():
    germs = [[1.0,1.0],[-1.0,2.0],[0.0,0.0],[2.0,1.0]]
    grid = adjust_grid(germs,[15,15])
    print grid
    plt.figure(2)
    for i in xrange(len(grid)):
        for j in xrange(len(grid[0])):
            if grid[i][j] == 0:
                plt.plot(float(i)*10.0/15.0-5.0,float(j)*10.0/15.0-5.0,'ro')
            elif grid[i][j] == 1:
                plt.plot(float(i)*10.0/15.0-5.0,float(j)*10.0/15.0-5.0,'go')
            elif grid[i][j] == 2:
                plt.plot(float(i)*10.0/15.0-5.0,float(j)*10.0/15.0-5.0,'yo')
            else:
                plt.plot(float(i)*10.0/15.0-5.0,float(j)*10.0/15.0-5.0,'bo')
    plt.show()

#displayregions()

























