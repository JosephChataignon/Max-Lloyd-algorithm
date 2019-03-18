
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
        return 0

def gaussian(x,y):
    return math.exp((-(x**2)-(y**2))/2)

# function studied
def f(x,y):
    if x**2+y**2 > 10:
        return 0
    return gaussian(x,y)

# distribution studied
def random_distrib():
    #return [random.uniform(-1,1),random.uniform(-1,1)]
    return [random.gauss(0,1),random.gauss(0,1)]

# squared distance between 2 points
def sqdistance(point1,point2):
    return (float(point1[0])-float(point2[0]))**2 + (float(point1[1])-float(point2[1]))**2


# computes mean squared error on R
def MSE(germs,boundaries,grid_dim):
    print "MSE"
    s = 0
    for i in xrange(grid_dim[0]):
        for j in xrange(grid_dim[1]):
            s = s + sqdistance( [i*10/grid_dim[0]-5,j*10/grid_dim[1]-5] , germs[grid[i][j]] ) / (grid_dim[0]*grid_dim[1])
    return s

# k is the index of the germ to adjust
def centroid(grid,k,grid_dim):
    print "centroid"
    Cx = integrate.quad(lambda z: Sx(z,boundaries), -10, 10)[0]
    if Cx != 0:
        Cx = integrate.quad(lambda z: z*Sx(z,boundaries), -10, 10)[0] / Cx
    else:
        Cx = 0
    Cy = integrate.quad(lambda z: Sy(z,boundaries), -10, 10)[0]
    if Cy != 0:
        Cy = integrate.quad(lambda z: z*Sy(z,boundaries), -10, 10)[0] / Cy
    else:
        Cy = 0
    return [Cx,Cy]

#germs is an array containing all the germs
def adjust_grid(germs,grid_dim):
    print "grid"
    grid = np.array([[0]*grid_dim[0]]*grid_dim[1])
    for i in xrange(len(grid)):
        for j in xrange(len(grid[0])):
            point = [float(i)*10.0/float(grid_dim[0])-5.0 , float(j)*10.0/float(grid_dim[1])-5.0]
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
def maxlloyd(germs,error_threshold,grid_dim):
    print "test fonction  maxlloyd"
    e = error_threshold+1
    error = [] # store the evolution of MSE through the loop
    c = 0 # counts the number of executions of the loop
    grid = [[0]*grid_dim[0]]*grid_dim[1]
    while e > error_threshold and c < 10:
        c = c+1
        print c
        if c%2 == 1:
            # adjust grid
            print "adjust grid"
            grid = adjust_grid(germs,grid_dim)
        else:
            # adjust germs
            print "adjust germs"
            for i in xrange(len(germs)):
                germs[i] = centroid(grid,i,grid_dim)
        e = MSE(germs,grid,grid_dim)
        error.append(e)
    return germs,error


# Test of maxlloyd function
def test_maxlloyd():
    g = [[1,1],[0,0],[2,1],[-1,2],[2,-1],[-4,-4]]
    g2,error = maxlloyd(g,0.1,[10,10])
    print g2,b2
    plt.plot(error)
    plt.show()
    return g2,b2

#test_maxlloyd()



# test function for visualising the regions of the grid after adjust_grid function
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



























