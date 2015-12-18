# hw07
  Dennis Scheglov   Joshua Pelton
  
  In this project, the goal was to use a recurrence relation to integrate (x^n)e^-x)dx from 0 to 1.
  The idea was that by performing a backwards recursion, one could avoid the emergent error in floating point values returned by
  the function from n = 0 to n = 100.
  The program also compares the performance of the recursive function with that of a high-quality general gsl integrator with 
  with machine E = 10^-9.
  Another goal was to compare the two functions in terms of the times they took per call, in a similar fashion to Midterm 1; by
  calling them multiple times and then dividing by the number of calls for greater accuracy.
  
## Results
   The general integrator took 1.1362176e-3 seconds per call.
   The recursive integration took 6.781900e-7 seconds per call.
   We found that the recursive integrator was 1.675454e03 times faster than the general integrator.
