#ifndef GEOM_LA_H_INCLUDED
#define GEOM_LA_H_INCLUDED

// Basic linear algebra functions for spatial geometry

unsigned geom_imin2f(const float  v[2]);
unsigned geom_imin2d(const double v[2]);
unsigned geom_imin3f(const float  v[3]);
unsigned geom_imin3d(const double v[3]);
unsigned geom_imin4f(const float  v[4]);
unsigned geom_imin4d(const double v[4]);
unsigned geom_imax2f(const float  v[2]);
unsigned geom_imax2d(const double v[2]);
unsigned geom_imax3f(const float  v[3]);
unsigned geom_imax3d(const double v[3]);
unsigned geom_imax4f(const float  v[4]);
unsigned geom_imax4d(const double v[4]);

// 2, 3, and 4D norms and normalization
float  geom_norm2f(const float  v[2]);
double geom_norm2d(const double v[2]);
float  geom_norm3f(const float  v[3]);
double geom_norm3d(const double v[3]);
float  geom_norm4f(const float  v[4]);
double geom_norm4d(const double v[4]);
float  geom_normalize2f(float  v[2]);
double geom_normalize2d(double v[2]);
float  geom_normalize3f(float  v[3]);
double geom_normalize3d(double v[3]);
float  geom_normalize4f(float  v[4]);
double geom_normalize4d(double v[4]);

float  geom_dot2f(const float  a[2], const float  b[2]);
double geom_dot2d(const double a[2], const double b[2]);
float  geom_dot3f(const float  a[3], const float  b[3]);
double geom_dot3d(const double a[3], const double b[3]);
float  geom_dot4f(const float  a[4], const float  b[4]);
double geom_dot4d(const double a[4], const double b[4]);

// Scalar cross product in 2D
float  geom_cross2f(const float  a[2], const float  b[2]);
double geom_cross2d(const double a[2], const double b[2]);

// Cross product in 3D
void geom_cross3f(const float  a[3], const float  b[3], float  result[3]);
void geom_cross3d(const double a[3], const double b[3], double result[3]);

// Given a vector (a), make two other vectors b and c such that
// a, b, c forms a right-handed orthogonal basis. b and c are normalized.
void geom_maketriad3f(const float  a[3], float  b[3], float  c[3]);
void geom_maketriad3d(const double a[3], double b[3], double c[3]);

// All matrices are stored in column major order:
//   [ m[0]  m[2] ]
//   [ m[1]  m[3] ]

// Overwrites v with the matrix-vector product m.v
void geom_matvec2f(const float  m[4],  const float  x[2], float  y[2]);
void geom_matvec2d(const double m[4],  const double x[2], double y[2]);
void geom_matvec3f(const float  m[9],  const float  x[3], float  y[3]);
void geom_matvec3d(const double m[9],  const double x[3], double y[3]);
void geom_matvec4f(const float  m[16], const float  x[4], float  y[4]);
void geom_matvec4d(const double m[16], const double x[4], double y[4]);

void geom_matTvec2f(const float  m[4],  const float  x[2], float  y[2]);
void geom_matTvec2d(const double m[4],  const double x[2], double y[2]);
void geom_matTvec3f(const float  m[9],  const float  x[3], float  y[3]);
void geom_matTvec3d(const double m[9],  const double x[3], double y[3]);
void geom_matTvec4f(const float  m[16], const float  x[4], float  y[4]);
void geom_matTvec4d(const double m[16], const double x[4], double y[4]);

// Overwrites b with the matrix-matrix product a.b
void geom_matmat2f(const float  a[4],  const float  b[4], float  c[4]);
void geom_matmat2d(const double a[4],  const double b[4], double c[4]);
void geom_matmat3f(const float  a[9],  const float  b[9], float  c[9]);
void geom_matmat3d(const double a[9],  const double b[9], double c[9]);
void geom_matmat4f(const float  a[16], const float  b[16], float  c[16]);
void geom_matmat4d(const double a[16], const double b[16], double c[16]);

// Inverts the matrix m in-place (assumes the matrix is well conditioned)
void geom_matinv2f(float  m[4]);
void geom_matinv2d(double m[4]);
void geom_matinv3f(float  m[9]);
void geom_matinv3d(double m[9]);
void geom_matinv4f(float  m[16]);
void geom_matinv4d(double m[16]);

// Computes the SVD of m = u.diag(s).vt
void geom_matsvd2d(const double m[4], double u[4], double s[2], double vt[4]);

// Solves the quadratic equation
//   a*x^2 + 2*b*x + c == 0
// Note the extra factor of 2 on the linear term.
int geom_quadraticd(
	const double a,
	const double b,
	const double c,
	double root[2]
);

#endif // GEOM_LA_H_INCLUDED
