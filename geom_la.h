#ifndef GEOM_LA_H_INCLUDED
#define GEOM_LA_H_INCLUDED

// Basic linear algebra functions for spatial geometry

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
void geom_matvec2f(const float  m[4],  float  v[2]);
void geom_matvec2d(const double m[4],  double v[2]);
void geom_matvec3f(const float  m[9],  float  v[3]);
void geom_matvec3d(const double m[9],  double v[3]);
void geom_matvec4f(const float  m[16], float  v[4]);
void geom_matvec4d(const double m[16], double v[4]);

// Overwrites b with the matrix-matrix product a.b
void geom_matmat2f(const float  a[4],  float  b[4]);
void geom_matmat2d(const double a[4],  double b[4]);
void geom_matmat3f(const float  a[9],  float  b[9]);
void geom_matmat3d(const double a[9],  double b[9]);
void geom_matmat4f(const float  a[16], float  b[16]);
void geom_matmat4d(const double a[16], double b[16]);

// Inverts the matrix m in-place (assumes the matrix is well conditioned)
void geom_matinv2f(float  m[4]);
void geom_matinv2d(double m[4]);
void geom_matinv3f(float  m[9]);
void geom_matinv3d(double m[9]);
void geom_matinv4f(float  m[16]);
void geom_matinv4d(double m[16]);

#endif // GEOM_LA_H_INCLUDED
