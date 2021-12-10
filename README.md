# Precise International Roughness Index (IRI) calculation
Matlab and Python codes for the computation of International Roughness Index (IRI) provided with our paper 
F. Sroubek, M. Sorel, J. Zak: Precise International Roughness Index Calculation. Int. J. Pavement Res. Technol. (2021).
If you use this code in your academic work, please cite our paper properly according
https://doi.org/10.1007/s42947-021-00097-z . The paper can be found also on our page 
http://www.utia.cas.cz/biblio?pub=0545847 . The Python version also contains a command line tool to work
directly on files with profile data.

The core function iri (in iri.m for Matlab and iri.py for Python) contains the implementation of three methods: original
Sayers' method, using numerical solver ode45 and the method
from our paper. Our method works directly for irregular sampling and regular sampling
with shifted starting position. In the Python version, we have not implemented 
the second option (numerical solver ode45), which is too slow to be used in practice. 

The python version can be also called as a command-line script with similar
interface. 
  
There are two example profiles in the repository:
- profile_1.txt ... regular sampling 0.25m 
- profile_2.txt ... irregular sampling of the same profile (artifically computed by interpolation from profile_1.txt)

You can run four Matlab examples
- example_1.m ... regular sampling, no overlap, 20m IRI segments
- example_2.m ... regular sampling, step 0.5m, 20m IRI segments 
- example_3.m ... irregular sampling, no overlap, 20m IRI segments
- example_4.m ... irregular sampling, step 0.5m, 20m IRI segments

To run the Python version, please make sure that you have installed the packages numpy, scipy and matplotlib
(in requirements.txt). The command line interface has the same parameters as the core function,
except that both the input and IRI output are contained in text/csv files. It is also possible to plot
an IRI graph to a separate graphical file. For the description of command line parameters, see its help
by running python iri.py -h 

Examples of commandline use:
    
    python iri.py ../profile_1.txt testdir/iri_file.txt -plot_file testdir/plot_file.pdf
    python iri.py ../profile_2.txt testdir/iri_file2.txt -plot_file testdir/plot_file2.svg -step 0.25


--------------------------------------------------------
Parameters of the core iri function are the same for both Matlab and Python
--------------------------------------------------------


function IRI = iri(Y, segment_length, start_pos, step, box_filter, method)

Input:

Y ... road profile, Nx2 vector [x,y], both x and y must be in meters.

segment_length ... length of the IRI segment (typically 20 or 100 meters)

start_pos ... starting position of IRI segments, if empty, set to the first sample

step ... if empty, false or 0, IRI is computed in segments without overlap,
              otherwise, the step is the shift of consecutive IRI segments (can be
              set for example to sampling step for regular sampling);
              overlap cannot be logical true.

box_filter ... true/false (impl. value true), if true,
              profile is first averaged with a box filter of length 0.25m 
              as recommended in Sayers' paper for sampling intervals
              shorter than 0.25m to better represent the way in which 
              the tire of a vehicle envelops the ground.

method ...    0 - Sayers' implementation, for regular sampling only,
              1 - numerical solver ode45 (not implemented in Python),
              2(default) - our semi-analytical solution (Sroubek & Sorel)

Output:

IRI ... IRI vector calculated in segments; matrix with rows
          [star_pos, end_pos, IRI, std_of_IRI], where std_of_IRI is 
          standard deviation of IRI within each segment, IRI is in mm/m [equivalent to m/km]


