
/* Return a positive value if the point pd lies inside the     */
/* circle passing through pa, pb, and pc; a negative value if  */
/* it lies outside; and zero if the four points are cocircular.*/
/* The points pa, pb, and pc must be in counterclockwise       */
/* order, or the sign of the result will be reversed.          */
double orient2d(const double *pa, const double *pb, const double *pc);

/* Return a positive value if the point pd lies below the      */
/* plane passing through pa, pb, and pc; "below" is defined so */
/* that pa, pb, and pc appear in counterclockwise order when   */
/* viewed from above the plane.  Returns a negative value if   */
/* pd lies above the plane.  Returns zero if the points are    */
/* coplanar.  The result is also a rough approximation of six  */
/* times the signed volume of the tetrahedron defined by the   */
/* four points.                                                */
double orient3d(const double *pa, const double *pb, const double *pc, const double *pd);

/* Return a positive value if the point pd lies inside the     */
/* circle passing through pa, pb, and pc; a negative value if  */
/* it lies outside; and zero if the four points are cocircular.*/
/* The points pa, pb, and pc must be in counterclockwise       */
/* order, or the sign of the result will be reversed.          */
double incircle(const double *pa, const double *pb, const double *pc, const double *pd);

/* Return a positive value if the point pe lies inside the     */
/* sphere passing through pa, pb, pc, and pd; a negative value */
/* if it lies outside; and zero if the five points are         */
/* cospherical.  The points pa, pb, pc, and pd must be ordered */
/* so that they have a positive orientation (as defined by     */
/* orient3d()), or the sign of the result will be reversed.    */
double insphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe);


/*

Robust Predicates on Pentium CPUs

Intel CPUs need to be configured to correctly execute my code for arbitrary
floating-point precision arithmetic and robust geometric predicates. The
problem with these CPUs in their default state is that they use extended
precision internal floating-point registers, which defeat Dekker's standard
techniques for computing the roundoff error of a floating-point operation.
These processors can be configured to round internally stored floating-point
values to single or double precision, but the commands for doing so are
compiler-dependent.

Depending on what C compiler and operating system you use with your PC, you
might need to figure out your compiler's procedures for adjusting the
floating-point control register, and add them to the file predicates.c or
triangle.c. The correct place to add them is at the beginning of the
procedure exactinit(). You will probably also need to add an #include
statement somewhere early in the file.

If you are compiling with gcc under LINUX, be sure the LINUX symbol is
defined during the compilation. If you are compiling with Microsoft C, be
sure the CPU86 symbol is defined during the compilation. These symbols cause
the inclusion of the following code for configuring the floating-point
registers, which I suspect works, but I'm not really sure.

For gcc running under Linux, the statements are as follows. The following
line appears with the other #include statements.

#include <fpu_control.h>
The following lines appear at the beginning of the procedure exactinit().
  int cword;

#ifdef SINGLE
  cword = 4210;                 // set FPU control word for single precision
#else // not SINGLE
  cword = 4722;                 // set FPU control word for double precision
#endif // not SINGLE
  _FPU_SETCW(cword);
For Microsoft C, the incantations are as follows.
#include <float.h>

#ifdef SINGLE
_control87(_PC_24, MCW_PC);     // set FPU control word for single precision
#else // not SINGLE
_control87(_PC_53, MCW_PC);     // set FPU control word for double precision
#endif // not SINGLE

If you have any corrections to this information, or know how to configure Intel
or other FPUs under other compilers or operating systems, please send me email at
[omitted].
Alternatively, the following might also work with gcc, even if you're not using
Linux. (The definition of set_ctrlword() must occur before the definition of
exactinit().)

void set_ctrlword(v)
int v;
{
  asm("fldcw %0" :: "m" (v));
}

#ifdef SINGLE
  set_ctrlword(4210);           // set FPU control word for single precision
#else // not SINGLE
  set_ctrlword(4722);           // set FPU control word for double precision
#endif // not SINGLE

*/

