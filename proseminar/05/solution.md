# Exercise 1

## Measure the execution time for various problem sizes. What can you observe?

![Alt text](task1/plot/VaryNumberBodies.png "Title")

![Alt text](task1/plot/VaryNumberSteps.png "Title")

For varying number of bodies, as well as varying number of steps we have a linear dependency between the problem size and the execution time. 

# Exercise 2

## What optimization methods can you come up with in order to improve the performance of Exercise 1?
- We could parallize the process (computation of bodies can be spilt up to multiple cores)

## What parallelization strategies would you consider for Exercise 1 and why?
- A node is responsible for n bodies to calculate and at each time step we synchronize the entire data (Very high communication overhead => not lucrative) 
- Split up the space into multiple cells. Each node calulates the body if its center is inside of its responsible field. (Still high communication overhead)
- Barnes-Hut-Algorithmus (https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation): The simulation volume is usually divided up into cubic cells via an octree (in a three-dimensional space), so that only particles from nearby cells need to be treated individually, and particles in distant cells can be treated as a single large particle centered at the cell's center of mass (or as a low-order multipole expansion). This can dramatically reduce the number of particle pair interactions that must be computed.

 
