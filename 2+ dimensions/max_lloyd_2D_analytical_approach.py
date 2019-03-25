# -*- coding: utf-8 -*-
#
# Max-Lloyd algorithm for finding an optimal quantizer
# in dimension 2
#
# Warning!
# This approach didn't work because of the calculations of integrals, the code
# is still not functionning


import math
import random
import scipy
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

# indicator function: returns 1 if x is in the boundaries, 0 otherwise
def is_part_of_region(x, y, boundaries):
    for b in boundaries:
        if b[0]*x + b[1]*y > b[2] :
            return 0
    return 1

# same as f, but the returned value is zero if (x,y) is out of the boundaries
# here the array boundaries only contains the boundaries for one region
def f_with_constraints(u,v,boundaries):
    if is_part_of_region(u,v, boundaries) == 1:
        return f(u,v)
    else:
        return 0

# returns the measure of intersection between the zone inside boundaries and the line x=z
def Sx(z,boundaries):
    return integrate.quad(lambda t: is_part_of_region(z,t,boundaries), -float('Inf'), float('Inf'))[0]

# returns the measure of intersection between the zone inside boundaries and the line y=z
def Sy(z,boundaries):
    return integrate.quad(lambda t: is_part_of_region(t,z,boundaries), -float('Inf'), float('Inf'))[0]

# computes MSE inside one region
# here the array boundaries only contains the boundaries for the region of x
def region_MSE(x,boundaries):
    y = lambda u,v: f_with_constraints(u,v,boundaries) * ((u-x[0])**2+(v-x[1])**2)
    return integrate.dblquad( y, -float('Inf'), float('Inf'), -float('Inf'), float('Inf'))[0]

# computes mean squared error on R
def MSE(germs,boundaries):
    print "MSE"
    s = 0
    for i in xrange(len(germs)):
        s = s + region_MSE(germs[i],boundaries[i])
    return s

# boundaries is an array containing the constraints for one region
def centroid(boundaries):
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

# germs is an array containing all the germs
def adjust_boundaries(germs,k):
    print "boundaries"
    bounds = []
    main_germ = germs[k]
    for i in xrange(len(germs)):
        if i != k:
            a1 = 2*( germs[i][0] - main_germ[0] )
            a2 = 2*( germs[i][1] - main_germ[1] )
            b  = germs[i][0]**2 + germs[i][1]**2 - main_germ[0] - main_germ[1]
            bounds.append([a1,a2,b]) 
    return bounds

# boundaries is an array containing the initial decision boundaries (array of n-1 constraints)
# constraints' model is [a1,a2,b] parameters of the equation a1*x + a2*y <= b
# germs is an array containing the germ for each region, on the model [x,y]
# error_threshold is the threshold to reach for the algorithm to stop
def maxlloyd(germs,boundaries,error_threshold):
    print "test fonction  maxlloyd"
    e = error_threshold+1
    error = [] # store the evolution of MSE through the loop
    c = 0 # counts the number of executions of the loop
    while c < 10:#e > error_threshold and c < 10:
        c = c+1
        print c
        if c%2 == 1:
            # adjust boundaries
            print "adjust boundaries"
            for i in xrange(len(boundaries)):
                boundaries[i] = adjust_boundaries(germs,i)
        else:
            # adjust germs
            print "adjust germs"
            for i in xrange(len(germs)):
                germs[i] = centroid(boundaries[i])
        e = MSE(germs,boundaries)
        error.append(e)
    return germs,boundaries,error


# Test for maxlloyd function
def test_maxlloyd():
    g = [[1,1],[0,0],[2,1],[-1,2],[2,-1],[-4,-4]]
    b = [[[0,0,0]]*5]*6
    g2,b2,error = maxlloyd(g,b,0.1)
    print g2,b2
    plt.plot(error)
    plt.show()
    return g2,b2

#test_maxlloyd()




def test_is_part_of_region():
    g = [[1.0,1.0],[0.0,0.0],[2.0,1.0],[-1.0,2.0],[2.0,-1.0],[-4.0,-4.0]]
    b = [[[0.0,0.0,0.0]]*5]*6
    for i in xrange(len(b)):
        b[i] = adjust_boundaries(g,i)
    tab = []
    
    for i in xrange(100):
        for j in xrange(100):
            z = is_part_of_region(i/10-5,j/10-5,b[5])
            tab.append([i/10-5,j/10-5,z])
    plt.figure(2)
    for point in tab:
        if point[2] == 0:
            plt.plot(point[0],point[1],'bo')
        else:
            plt.plot(point[0],point[1],'ro')
    plt.show()

































