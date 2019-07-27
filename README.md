# Max-Lloyd-algorithm
This repository contains an implementation of Lloyd-Max algorithm in Python. It is divided as follows:
- the folder "1 dimension" contains an algorithm for 1-dimensional problems, as well as a verification table
 of values to check the accuracy of results. It is the proper Max-Lloyd algorithm.
- the folder "2+ dimensions" contains 3 code files corresponding to 3 different approaches of the problem in
 2 dimensions or more, with some of the results obtained with the mixed approach. It should in fact be called 
 Linde-Buzo-Gray algorithm as Max-Lloyd was initially only for one dimension.

I started  with the "analytical approach", which was suposed to give the result only with analytical calculations (and only in 2 dimensions).
Hovewer it met an obstacle when infinite integrals calculations were needed. Thus, I switched to a completely numerical approach,
where space  is modelised by a grid of points, that worked (in 2 dimensions) but with bad efficiency.
Finally I created a mix of both approaches by adding some numerical calculations on an analytical modelisation, and adapting 
this code for higher dimensions.

Therefore if you are looking for a functional implementation, you should use the "mixed approach" code. It comes with several testing functions included
 in the second half of the file.
