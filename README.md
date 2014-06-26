# CFD
## A fluid dynamics & radiation transport simulation code

******************************


**Overview**
A 3D Eulerian fixed grid fluid dynamics code. The grid is a collection of finite elements, called grid cells, that 
each hold fluid state information. The hydrodynamics are solved using a rotated hybrid HLLC-HLL Riemann solver
to calculate fluxes on each grid cell face. Ionization from point source radiation is implicitly solved and the 
column densities required for this are calculated via the method of short characteristics. Heating/cooling from 
atomic processes is calculated using the approximate functions in Henney (2009).

* **Example Usage**
TODO.

* **Getting Started**
Requires Eigen: "a template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms"

http://eigen.tuxfamily.org

and TinyXml: "a simple, small, C++ XML parser"

http://www.grinninglizard.com/tinyxml/

* **Advanced Usage**
TODO.

* **Goals**
None at the moment.

* **Developer info**
TODO.

* **Colophon**
TODO: references.

## Supporting Documentation
TODO.
