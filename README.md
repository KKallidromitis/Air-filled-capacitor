# Air-filled-capacitor
The matlab code is used to estimate the capacitance C per unit length of a long air-filled capacitor based on a mesh size (h) which solves a set of differential equations.Alfa is the over relaxation parameter which determines how fast the Finite Differences (FD) method converges to a value. Each mesh size h has only one value for alpha which represents the minimum number of iterations for the solution. 

Estimating_optimal_a.m initializes alpha to be close to 2 then runs 10000 times from (1.9999 to 1.0001) calculating each iteration (iteration loop inside while loop) for every alpha and placing them inside the matrices x1, y1. Then plots a graph with all values of alpha (1<α <2) and their corresponding iterations for h=0.5.

test.m is used to solve the differential equations and compare the solutions with the experimental result from the gradient descent algorithm.
