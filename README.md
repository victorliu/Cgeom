Cgeom
=====
A set of computational geometry routines in C intended for engineering simulation.

* `geom_la`: Basic linear algebra routines for 2/3/4 dimensions.
* `geom_poly`: Polygon and polyhedron routines in 2d and 3d. Depends on `geom_la`.
* `geom_shapes`: Basic shapes, bounds, and predicates. Depends on `geom_la` and `geom_poly`.
* `geom_bvh`: Bounding volume hierarchy object.
* `geom_shapeset`: High level operations on a collection of shapes. Depends on `geom_shapes` and `geom_bvh`.
* `geom_predicates`: Robust geometric predicates.
* `geom_circum`: Circumcenter routines. Depends on `geom_predicates` for high precision.
