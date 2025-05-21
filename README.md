# N body simulation project for CSC_3S005_EP

## Prerequisites

Input data: The description of $N$ bodies $b_1,\ldots,b_N$: mass $m_i$, initial coordinates $x_i,y_i$, initial velocity $u_i,v_i$

The force two bodies $b_i, b_j$ excert on each other is given by the formula:
$$F_{ij} = \frac{G m_i m_j}{(x_i-x_j)^2+(y_i-y_j)^2}$$
where G is the gravitational constant.

## Distributing the work
we have 3 main thing to work on:
*Main thread*: create threads and make sure they are sincronysed in the manner we want (every step or every T steps), save the results

*Calculation threads*: Take a chunck of a matrix and take a matrix of distances and then they are going to calculate on that chunck of the matrix the updates to the velocities. They take the timestep at which they start and then we synchronise the threads that are to far behind.

*(?)Body calculations*: Work on the new velocities and distances of each body. They need to be time step synchronised.

We need to use locks.

Parameters: time $t$, number of bodies $N$ (position $x_i, y_i$, velocity, mass), $\varepsilon$ timestamp sync, $th$ number of threads

Basic structure: Gabriel
Main thread: Daniela
Force compitation threads: Aya (Gabriel as backup)
Body computation threads: Gabriel (Daniela as backup)
Graphics: Aya (Daniela as backup)
