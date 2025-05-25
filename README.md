# N body simulation project for CSC_3S005_EP

## Prerequisites

Input data: 
- The description of $N$ bodies $b_1,\ldots,b_N$: mass $m_i$ (in kilograms, kg), initial coordinates $x_i,y_i$ (in meters, m), initial velocity $u_i,v_i$ (in meters per second, m*s^-1). Apart from the number of bodies, each data point is optional and will be randomly initialised if not given (or set to 0 if only velocity is not given).
- The time (in seconds, s) each time step (or frame) takes. Longer times will result in much faster calculations but likely more innacurate ones.
- The time step synchronization $\varepsilon$ (an integer) which determines how often forces are *required* to be updated (by locking required threads until forces have been updated). Ensures the result is not too innacurate. A time step synchronization of 1 will require the forces to be recalculated every time step.
- The number of threads.

The force two bodies $b_i, b_j$ excert on each other is given by the formula:
$$F_{ij} = \frac{G m_i m_j}{(x_i-x_j)^2+(y_i-y_j)^2}$$
where $G$ is the gravitational constant (equal to 6.6743 × 10^-11).

Two global matrices:
- The first, the "Force" matrix, is a triangular/symmetrical matrix which stores at position $i, j$ the force (in Newtons, N, kg\*m\*s^-2) body $i$ exerts on body $j$. Position $i, j$ is further subdivided into two parts, for the force in the $x$ direction and the force in the $y$ direction.
- The second, the "Distance" matrix, is a triangular/symmetrical matrix which stores at position $i, j$ the distance (in meters, m) between body $i$ and body $j$.

## Distributing the work
There are 3 main parts:
- *Main Thread*: Create threads and make sure they are synchronised in the chosen way (every $\varepsilon$ steps), and save the results.
- *Force Calculation Threads*: Assign from the Force matrix a small chunck to each thread. Each thread then calculates using the Distance matrix the forces required for its chunck. Each thread reads the current time step $t$ when starting its update, and when done, stores that value as its latest "updated" time step before restarting. Calculation threads may need to grab locks if they are more than $\varepsilon$ time steps behind the current time step (although there is probably a better way to do this).
- *Body Calculation Threads*: Each thread is assigned a certain number of bodies. For each body, calculate its new velocity by summing the forces on the corresponding row of the Force matrix, then update the body's position and velocity. Then update the Distance matrix (optionally, depending on thread speeds, the Force Calculation Threads may update the Distance Matrix before calculating Forces, instead of Body Calculation Threads updating the Distance matrix). Body calculation threads must be time step synchronized (using locks if necessary).

Basic structure: Gabriel
Main thread: Daniela
Force computation threads: Aya (Gabriel as backup)
Body computation threads: Gabriel (Daniela as backup)
Graphics: Aya (Daniela as backup)
