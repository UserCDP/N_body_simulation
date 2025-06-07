# N body simulation project for CSC_3S005_EP

## Prerequisites

Input data:

- The description of $N$ bodies $b_1,\ldots,b_N$: mass $m_i$ (in kilograms, kg), initial coordinates $x_i,y_i$ (in meters, m), initial velocity $u_i,v_i$ (in meters per second, m\*s^-1). Apart from the number of bodies, each data point is optional and will be randomly initialised if not given (or set to 0 if only velocity is not given).
- The time (in seconds, s) each time step (or frame) takes. Longer times will result in much faster calculations but likely more innacurate ones.
- The number of time_steps to simulate.
- The number of threads.

The force two bodies $b_i, b_j$ excert on each other is given by the formula:
$$F_{ij} = \frac{G m_i m_j}{(x_i-x_j)^2+(y_i-y_j)^2}$$
where $G$ is the gravitational constant (equal to 6.6743 × 10^-11).

Two matrices on the heap:

- The "Force" matrix, is a triangular/symmetrical matrix which stores at position $i, j$ the force (in Newtons, N, kg\*m\*s^-2) body $i$ exerts on body $j$. Position $i, j$ is further subdivided into two parts, for the force in the $x$ direction and the force in the $y$ direction.
- The "All_Positions" matrix, is a matrix which stores at position $i, j$ the position of body $j$ at time step $i$, which is further divided into $x$ coordinate and $y$ coordinate.
The "Body" instances are also stored on the heap.

## Distributing the work

There are 3 main parts:

- _Main Thread_: Create threads and synchronize threads as well as save the results.
- _Force Calculation Threads_: Assign from the Force matrix a small chunck to each thread. Each thread then calculates the forces required for its chunck.
- _Body Calculation Threads_: Each thread is assigned a certain number of bodies. For each body, calculate its new velocity by summing the forces on the corresponding row of the Force matrix, then update the body's position and velocity.

Main thread and Basic Structure: Daniela (Gabriel helping)

Sequential and concurrent simulation code: Gabriel (Aya helping)

Graphics/Visualization: Aya (Daniela as backup)
