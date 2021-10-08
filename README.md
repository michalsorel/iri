# Precise International Roughness Index (IRI) calculation
Matlab code for the computation of International Roughness Index (IRI) provided with our paper 
F. Sroubek, M. Sorel, J. Zak: Precise International Roughness Index Calculation. Int. J. Pavement Res. Technol. (2021).
If you use this code in your academic work, please cite our paper properly according
https://doi.org/10.1007/s42947-021-00097-z . The paper can be found also on our page 
http://www.utia.cas.cz/biblio?pub=0545847 . We plan to add a Python version in near future, stay tuned.

There are two example profiles in the repository:
- profile_1.txt ... regular sampling 0.25m 
- profile_2.txt ... irregular sampling of the same profile 

You can run four examples
- example_1.m ... regular sampling, no overlap, 20m IRI segments
- example_2.m ... regular sampling, overlap 0.5m, 20m IRI segments 
- example_3.m ... irregular sampling, no overlap, 20m IRI segments
- example_4.m ... irregular sampling, overlap 0.5m, 20m IRI segments

The core function iri contains the implementation of three methods: original
Sayers method, using numerical solver ode45 and the method
from our paper. Our method works also for irregular sampling and regular sampling
with shifted starting position.


--------------------------------------------------------


function IRI = iri(Y, segment_length, start_pos, overlap, box_filter, method)

Input:

Y ... road profile, Nx2 vector [x,y], both x and y must be in meters.

segment_length ... length of the IRI segment (typically 20 or 100 meters)

start_pos ... starting position of IRI segments, if empty, set to the first sample

overlap ... if empty, false or 0, IRI is computed in segments without overlap,
              otherwise, the overlap is overlap of IRI segments (can be
              set for example to sampling step for regular profile sampling);
              overlap cannot be logical true.

box_filter ... true/false (impl. value true), if true,
              profile is first averaged with a box filter of length 0.25m 
              as recommended in Sayers' paper for sampling intervals
              shorter than 0.25m to better represent the way in which 
              the tire of a vehicle envelops the ground.

method ...    0 - Sayers' implementation, for regular sampling only,
              1 - numerical solver ode45,
              2(default) - our semi-analytical solution (Sroubek & Sorel)

Output:

IRI ... IRI vector calculated in segments; matrix with rows
          [star_pos, end_pos, IRI, std_of_IRI], where std_of_IRI is 
          standard deviation of IRI within each segment, IRI is in mm/m [equivalent to m/km]


