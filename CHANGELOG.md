# [0.2]
* Changed Boundary class from one that describes individual GridJoin-like objects (that connect to the GridCells on the boundary) to one that contains all the cells that connect onto a face of the main grid. Each of these GridCell objects (or ghostcells) are linked to other GridCells, the number of which depends on the spatial and temporal order of integration.
* The member function Boundary::applyBC() now loops over all GridCells in the Boundary and applies the boundary condition to them.
