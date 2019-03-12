
# Max-Lloyd algorithm for finding the optimal quantizer
# in dimension 2

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
    return math.exp(-x**2/2-y**2/2)

# function studied
def f(x,y):
    return gaussian(x,y)

# distribution studied
def random_distrib():
    #return [random.uniform(-1,1),random.uniform(-1,1)]
    return [random.gauss(0,1),random.gauss(0,1)]

# indicator function: returns 1 if x is in the boundaries, 0 otherwise
def is_part_of_region(x, y, boundaries):
    for b in boundaries:
        print b
        if b[0]*x + b[1]*y > b[2] :
            return 0
    return 1

# same as f, but the returned value is zero if (x,y) is out of the boundaries
# here the array boundaries only contains the boundaries for the concerned region
def f_with_constraints(x,y,boundaries):
    if is_part_of_region(x, y, boundaries) == 1:
        return f(x,y)
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
    s = 0
    for i in xrange(len(germs)):
        s = s + region_MSE(germs[i],boundaries[i])
    return s

# boundaries is an array containing the constraints for one region
def centroid(boundaries):
    Cx = integrate.quad(lambda z: z*Sx(z,boundaries), -float('Inf'), float('Inf'))[0]
    Cx = Cx / integrate.quad(lambda z: Sx(z,boundaries), -float('Inf'), float('Inf'))[0]
    Cy = integrate.quad(lambda z: z*Sy(z,boundaries), -float('Inf'), float('Inf'))[0]
    Cy = Cy / integrate.quad(lambda z: Sy(z,boundaries), -float('Inf'), float('Inf'))[0]
    return [Cx,Cy]


def adjust_boundaries(germs,k):
    bounds = []
    main_germ = germs[k]
    for i in xrange(len(germs)):
        a1 = 2*( germs[i][0] - main_germ[0] )
        a2 = 2*( germs[i][1] - main_germ[1] )
        b  = germs[i][0]**2 + germs[i][1]**2 - main_germ[0] - main_germ[1]
        bounds.append([a1,a2,b]) 
    return bounds

# boundaries is an array containing the initial decision boundaries (array of n-1 constraints)
# constraints' model is [a1,a2,b] parameters of the equation a1*x + a2*y <= b
# germs is an array containing the germ for each region, on the model [x,y]
# error_threshold is the threshold to reach for the algorithm to stop
def maxlloyd(boundaries,germs,error_threshold):
    e = MSE(boundaries,germs)
    error = [e] # store the evolution of MSE through the loop
    c = 0 # counts the number of executions of the loop
    while e > error_threshold and c < 300:
        c = c+1
        if c%2 == 1:
            # adjust boundaries
            for i in xrange(len(boundaries)):
                boundaries[i] = adjust_boundaries(germs,i)
        else:
            # adjust germs
            for i in xrange(len(germs)):
                germs[i] = centroid(boundaries[i])
        e = MSE(boundaries,germs)
        error.append(e)
    return germs,boundaries,error

[0,0,0],
# Test of maxlloyd function
def test_maxlloyd():
    g = [[1,1],[0,0],[2,1],[-1,2],[2,-1],[-4,-4]]
    b = [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    g2,b2,error = maxlloyd(g,b,0.1)
    print x2,t2
    plt.plot(error)
    plt.show()
    return x2,t2

test_maxlloyd()

def estimate(x,t,value):
    for i in xrange(len(t)):
        if t[i] > value:
            return x[i]
    return x[-1]

# Plot of average error
def plot_avg_error(N):
    x,t = test_maxlloyd()
    avg_E = []
    realizations = []
    square_error = []
    for i in xrange(N):
        realizations.append(random_distrib())
        square_error.append((realizations[-1] - estimate(x,t,realizations[-1]))**2)
        avg_E.append(sum(square_error)/len(square_error))
    plt.figure(2)
    plt.plot(avg_E)
    plt.show()

#plot_avg_error(20000)

































