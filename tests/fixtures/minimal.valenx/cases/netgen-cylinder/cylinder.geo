## Tiny Netgen CSG fixture — a unit-radius cylinder of length 2 along Z.
##
## The Netgen adapter copies this file into the run workdir before
## invoking `netgen -batchmode -geofile=cylinder.geo`. Kept small
## enough to round-trip through git without LFS.

algebraic3d

solid cyl = cylinder (0, 0, -1; 0, 0, 1; 0.5)
       and plane (0, 0, -1; 0, 0, -1)
       and plane (0, 0,  1; 0, 0,  1);

tlo cyl;
