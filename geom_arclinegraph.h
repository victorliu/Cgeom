/* Given a simply connected polyline (not closed) defined by a
 * sequence of n points (and so there are n-1 segments), produce
 * a closed, positively (CCW) oriented polygon that is the offset
 * curve of the polyline with a specified total thickness.
 * This function assumes the input does not have adjacent segments
 * that double back on themselves.
 */
int polyline_thicken(unsigned int n, const double *vin, double thickness, double *vout);

/* Given a planar graph with edges that are either line segments or
 * arcs, produce an offset polygon around the graph with a specified
 * total thickness. The function does not check for self-intersections
 * in either input or output. The input segments/arcs are specified
 * by their endpoint indices in a list of vertex positions v, and the
 * arc bulge factors are in g.
 *
 * Join type:
 *   0: Simple miter     1: Bevel           2: Rounded
 *  
 *      \  \  \           \  \  \             \  \  \
 *    ---*  \  \        ---*  \  \          ---*  \  \
 *    -------*  \       -------*  /         -------*  ;
 *    -----------*      ---------/          ---------'
 *
 * The return arrays w may hold up to 6*n vertices, icd may hold up to
 * 6*n segments. The reason for this accounting is that beveled and
 * rounded corners produce an extra segment shared between adjacent
 * original segments, and ends of segments when capped with rounded
 * caps require two arcs to form. So the absolute worst case scenario
 * is a bent pair of segments with rounded corners, producing 9 output
 * segments. Otherwise, asymptotically, another worst case is a set of
 * completely disjoint (unconnected) segments with rounded caps
 * produces 6 output segments each (4 each for miter and bevel) and 6
 * output vertices each (4 each for miter and bevel).
 */
int arclinegraph_thicken(
	int join_type, const double *v,
	/* input arcs: */
	unsigned int n, const unsigned int *iab, const double *g,
	double thickness,
	/* outputs: */
	double *w, unsigned int *icd, double *h,
	int *nw, int *nseg
);
