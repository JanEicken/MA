Trajectory Computing Library
==========

Welcome to the initial release of the Trajectory Computing Library 

	libtrajcomp

We intend to bring together clean, template-driven C++ implementation of 
important algorithms in the field of trajectory computing in order to
facilitate research in this domain. 
The choice for template-driven C++ is due to the importance of fast computing
and memory effectiveness unreachable with, for example, Java or pure python and that
we can easily integrate it into important application environments including

*	R
*	GNU Octave
*	Python
*	Android / iOS

If you are interested in Java implementations, you can find some on our
Java trajectory computing sub-page, mainly contributed by students of mine:

http://www.trajectorycomputing.com/java

If you can give useful comments or even contribute algorithms or patches,
send them to me via email (for the first time, we might turn to a platform
for that, once a community grows...)

Table of Contents 
================

Tools
-----------
* [makestring: A flexible vector to string converter](tools/make_string.html)
* [progress: A simple progress bar for Linux terminals](tools/progress.html)
* [tictoc: A simple timer object](tools/tictoc.html)
* [matrix_resize: Using nested vectors to store matrices](tools/matrix_resize.html)



Basic Geometry Implementations
-----------------
* [default\_element\_distance: Euclidean Distance between Points](tools/default_element_distance.html)
* [default\_segment\_distance: Euclidean Distance between a Point and a Line Segment](tools/default_element_segment_distance.html)
* [wgs84\_element\_distance: WGS84 Distance between Points](tools/wgs84_element_distance.html)
* [wgs84\_segment\_distance: WGS84 Distance between a Point and a Line Segment](tools/wgs84_element_segment_distance.html)



Trajectory Computing Implementation
--------------------
* [trajectory: The main class for storing trajectories](trajectory.html)
Preoprocessing
---------------------
* [uniform\_select: Elementary Downsampling of Trajectories](preprocessing/uniform_select.html)




We like to thank all contributors to our dependencies for their hard work and hope, this library is helpful to you.