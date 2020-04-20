# iri
Code for the computation of International Roughness Index (IRI)

function [IRI] = iri_v2(Y, METHOD, VERBOSE, segmPar, avgPar)


Input:

Y ... road profile, Nx2 vector [x,y], Both x and y must be in meters.

METHOD ... 0 - Sayers' implementation, 1 - numerical solver ode45,
               2(default) - our semi-analytical solution

VERBOSE ... 0(default) - no plot; 1 - plot IRI graph 

segmPar ... parameters of segments; vector containing 3 parameters: 
               [segment_length(m), overlapping_segments(0/1), starting_position(m)] 
             default: [20,0,0]
			 
avgPar ... averaging parameter for smoothing profile; convolution with
           Gaussian function of variance proportional to 1/avgPar
            default: None 

Output:

IRI ... IRI vector calculated in segments; matrix with rows
           [star_pos, end_pos, IRI, std_of_IRI]
		IRI is in mm/m [equivalent to m/km]


