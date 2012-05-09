#include <math.h>
#include <string.h>

float geom_norm2f(const float v[2]){
	float x = fabsf(v[0]);
	float y = fabsf(v[1]);
	if(x < y){
		if(0 == x){
			return y;
		}
		x /= y;
		return y * sqrtf(1 + x*x);
	}else{
		y /= x;
		return x * sqrtf(1 + y*y);
	}
}
double geom_norm2d(const double v[2]){
	double x = fabs(v[0]);
	double y = fabs(v[1]);
	if(x < y){
		if(0 == x){
			return y;
		}
		x /= y;
		return y * sqrt(1 + x*x);
	}else{
		y /= x;
		return x * sqrt(1 + y*y);
	}
}

float geom_norm3f(const float v[3]){
	float a[3] = {fabsf(v[0]),fabsf(v[1]),fabsf(v[2])};
	float w = a[0];
	if(a[1] > w){ w = a[1]; }
	if(a[2] > w){ w = a[2]; }
	if(0 == w){
		return a[0] + a[1] + a[2];
	}else{
		a[0] /= w; a[1] /= w; a[2] /= w;
		w *= sqrtf(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
		return w;
	}
}

double geom_norm3d(const double v[3]){
	double a[3] = {fabs(v[0]),fabs(v[1]),fabs(v[2])};
	double w = a[0];
	if(a[1] > w){ w = a[1]; }
	if(a[2] > w){ w = a[2]; }
	if(0 == w){
		return a[0] + a[1] + a[2];
	}else{
		a[0] /= w; a[1] /= w; a[2] /= w;
		w *= sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
		return w;
	}
}

float geom_norm4f(const float v[4]){
	float a[4] = {fabsf(v[0]),fabsf(v[1]),fabsf(v[2]),fabsf(v[3])};
	float w = a[0];
	if(a[1] > w){ w = a[1]; }
	if(a[2] > w){ w = a[2]; }
	if(a[3] > w){ w = a[3]; }
	if(0 == w){
		return a[0] + a[1] + a[2] + a[3];
	}else{
		a[0] /= w; a[1] /= w; a[2] /= w; a[3] /= w;
		w *= sqrtf(a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]);
		return w;
	}
}

double geom_norm4d(const double v[4]){
	double a[4] = {fabs(v[0]),fabs(v[1]),fabs(v[2]),fabs(v[3])};
	double w = a[0];
	if(a[1] > w){ w = a[1]; }
	if(a[2] > w){ w = a[2]; }
	if(a[3] > w){ w = a[3]; }
	if(0 == w){
		return a[0] + a[1] + a[2] + a[3];
	}else{
		a[0] /= w; a[1] /= w; a[2] /= w; a[3] /= w;
		w *= sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3]);
		return w;
	}
}

float  geom_normalize2f(float  v[2]){
	float n = geom_norm2f(v);
	v[0] /= n; v[1] /= n;
	return n;
}
double geom_normalize2d(double v[2]){
	double n = geom_norm2d(v);
	v[0] /= n; v[1] /= n;
	return n;
}
float  geom_normalize3f(float  v[3]){
	float n = geom_norm3f(v);
	v[0] /= n; v[1] /= n; v[2] /= n;
	return n;
}
double geom_normalize3d(double v[3]){
	double n = geom_norm3d(v);
	v[0] /= n; v[1] /= n; v[2] /= n;
	return n;
}
float  geom_normalize4f(float  v[4]){
	float n = geom_norm4f(v);
	v[0] /= n; v[1] /= n; v[2] /= n; v[3] /= n;
	return n;
}
double geom_normalize4d(double v[4]){
	double n = geom_norm4d(v);
	v[0] /= n; v[1] /= n; v[2] /= n; v[3] /= n;
	return n;
}

void geom_cross3f(const float  a[3], const float  b[3], float  result[3]){
	result[0] = a[1]*b[2] - a[2]*b[1];
	result[1] = a[2]*b[0] - a[0]*b[2];
	result[2] = a[0]*b[1] - a[1]*b[0];
}
void geom_cross3d(const double a[3], const double b[3], double result[3]){
	result[0] = a[1]*b[2] - a[2]*b[1];
	result[1] = a[2]*b[0] - a[0]*b[2];
	result[2] = a[0]*b[1] - a[1]*b[0];
}

void geom_maketriad3f(const float  a[3], float  b[3], float  c[3]){
	float alen = geom_norm3f(a);
	const float an[3] = {
		a[0] / alen,
		a[1] / alen,
		a[2] / alen
	};
	if(fabsf(a[0]) > fabsf(a[1])){
		float invLen = 1. / geom_norm2f(an[0], an[2]);
		b[0] = -an[2] * invLen;
		b[1] = 0;
		b[2] = an[0] * invLen;
	}else{
		float invLen = 1. / geom_norm2f(an[1], an[2]);
		b[0] = 0;
		b[1] = an[2] * invLen;
		b[2] = -an[1] * invLen;
	}
	geom_cross3f(an, b, c);
}
void geom_maketriad3d(const double a[3], double b[3], double c[3]){
	double alen = geom_norm3d(a);
	const double an[3] = {
		a[0] / alen,
		a[1] / alen,
		a[2] / alen
	};
	if(fabs(a[0]) > fabs(a[1])){
		double invLen = 1. / geom_norm2d(an[0], an[2]);
		b[0] = -an[2] * invLen;
		b[1] = 0;
		b[2] = an[0] * invLen;
	}else{
		double invLen = 1. / geom_norm2d(an[1], an[2]);
		b[0] = 0;
		b[1] = an[2] * invLen;
		b[2] = -an[1] * invLen;
	}
	geom_cross3d(an, b, c);
}

void geom_matvec2f(const float  m[4],  float  v[2]){
	float v0 = v[0];
	v[0] = m[0] * v[0] + m[2] * v[1];
	v[1] = m[1] * v0   + m[3] * v[1];
}
void geom_matvec2d(const double m[4],  double v[2]){
	double v0 = v[0];
	v[0] = m[0] * v[0] + m[2] * v[1];
	v[1] = m[1] * v0   + m[3] * v[1];
}
void geom_matvec3f(const float  m[9],  float  v[3]){
	float u[3] = { v[0], v[1], v[2] };
	v[0] = m[0] * u[0] + m[3] * u[1] + m[6] * u[2];
	v[1] = m[1] * u[0] + m[4] * u[1] + m[7] * u[2];
	v[2] = m[2] * u[0] + m[5] * u[1] + m[8] * u[2];
}
void geom_matvec3d(const double m[9],  double v[3]){
	double u[3] = { v[0], v[1], v[2] };
	v[0] = m[0] * u[0] + m[3] * u[1] + m[6] * u[2];
	v[1] = m[1] * u[0] + m[4] * u[1] + m[7] * u[2];
	v[2] = m[2] * u[0] + m[5] * u[1] + m[8] * u[2];
}
void geom_matvec4f(const float  m[16], float  v[4]){
	float u[4] = { v[0], v[1], v[2], v[3] };
	v[0] = m[0] * u[0] + m[4] * u[1] + m[ 8] * u[2] + m[12] * u[3];
	v[1] = m[1] * u[0] + m[5] * u[1] + m[ 9] * u[2] + m[13] * u[3];
	v[2] = m[2] * u[0] + m[6] * u[1] + m[10] * u[2] + m[14] * u[3];
	v[3] = m[3] * u[0] + m[7] * u[1] + m[11] * u[2] + m[15] * u[3];
}
void geom_matvec4d(const double m[16], double v[4]){
	double u[4] = { v[0], v[1], v[2], v[3] };
	v[0] = m[0] * u[0] + m[4] * u[1] + m[ 8] * u[2] + m[12] * u[3];
	v[1] = m[1] * u[0] + m[5] * u[1] + m[ 9] * u[2] + m[13] * u[3];
	v[2] = m[2] * u[0] + m[6] * u[1] + m[10] * u[2] + m[14] * u[3];
	v[3] = m[3] * u[0] + m[7] * u[1] + m[11] * u[2] + m[15] * u[3];
}

void geom_matmat2f(const float  a[4],  float  c[4]){
	float b[4] = { c[0], c[1], c[2], c[3] };
	c[0] = a[0]*b[0] + a[2]*b[1];
	c[1] = a[1]*b[0] + a[3]*b[1];
	c[2] = a[0]*b[2] + a[2]*b[3];
	c[3] = a[1]*b[2] + a[3]*b[3];
}
void geom_matmat2d(const double a[4],  double c[4]){
	double b[4] = { c[0], c[1], c[2], c[3] };
	c[0] = a[0]*b[0] + a[2]*b[1];
	c[1] = a[1]*b[0] + a[3]*b[1];
	c[2] = a[0]*b[2] + a[2]*b[3];
	c[3] = a[1]*b[2] + a[3]*b[3];
}
void geom_matmat3f(const float  a[9],  float  c[9]){
	float b[9]; memcpy(b, c, 9*sizeof(float));
	c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
	c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
	c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
	c[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
	c[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
	c[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
	c[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
	c[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
	c[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
}
void geom_matmat3d(const double a[9],  double c[9]){
	double b[9]; memcpy(b, c, 9*sizeof(double));
	c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
	c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
	c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
	c[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
	c[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
	c[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
	c[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
	c[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
	c[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
}
void geom_matat4f(const float  a[16], float  c[16]){
	float b[16]; memcpy(b, c, 16*sizeof(float));
	c[0] = a[0]*b[0] + a[4]*b[1] + a[ 8]*b[2] + a[12]*b[3];
	c[1] = a[1]*b[0] + a[5]*b[1] + a[ 9]*b[2] + a[13]*b[3];
	c[2] = a[2]*b[0] + a[6]*b[1] + a[10]*b[2] + a[14]*b[3];
	c[3] = a[3]*b[0] + a[7]*b[1] + a[11]*b[2] + a[15]*b[3];
	c[4] = a[0]*b[4] + a[4]*b[5] + a[ 8]*b[6] + a[12]*b[7];
	c[5] = a[1]*b[4] + a[5]*b[5] + a[ 9]*b[6] + a[13]*b[7];
	c[6] = a[2]*b[4] + a[6]*b[5] + a[10]*b[6] + a[14]*b[7];
	c[7] = a[3]*b[4] + a[7]*b[5] + a[11]*b[6] + a[15]*b[7];
	c[ 8] = a[0]*b[8] + a[4]*b[9] + a[ 8]*b[10] + a[12]*b[11];
	c[ 9] = a[1]*b[8] + a[5]*b[9] + a[ 9]*b[10] + a[13]*b[11];
	c[10] = a[2]*b[8] + a[6]*b[9] + a[10]*b[10] + a[14]*b[11];
	c[11] = a[3]*b[8] + a[7]*b[9] + a[11]*b[10] + a[15]*b[11];
	c[12] = a[0]*b[12] + a[4]*b[13] + a[ 8]*b[14] + a[12]*b[15];
	c[13] = a[1]*b[12] + a[5]*b[13] + a[ 9]*b[14] + a[13]*b[15];
	c[14] = a[2]*b[12] + a[6]*b[13] + a[10]*b[14] + a[14]*b[15];
	c[15] = a[3]*b[12] + a[7]*b[13] + a[11]*b[14] + a[15]*b[15];
}
void geom_matmat4d(const double a[16], double c[16]){
	double b[16]; memcpy(b, c, 16*sizeof(double));
	c[0] = a[0]*b[0] + a[4]*b[1] + a[ 8]*b[2] + a[12]*b[3];
	c[1] = a[1]*b[0] + a[5]*b[1] + a[ 9]*b[2] + a[13]*b[3];
	c[2] = a[2]*b[0] + a[6]*b[1] + a[10]*b[2] + a[14]*b[3];
	c[3] = a[3]*b[0] + a[7]*b[1] + a[11]*b[2] + a[15]*b[3];
	c[4] = a[0]*b[4] + a[4]*b[5] + a[ 8]*b[6] + a[12]*b[7];
	c[5] = a[1]*b[4] + a[5]*b[5] + a[ 9]*b[6] + a[13]*b[7];
	c[6] = a[2]*b[4] + a[6]*b[5] + a[10]*b[6] + a[14]*b[7];
	c[7] = a[3]*b[4] + a[7]*b[5] + a[11]*b[6] + a[15]*b[7];
	c[ 8] = a[0]*b[8] + a[4]*b[9] + a[ 8]*b[10] + a[12]*b[11];
	c[ 9] = a[1]*b[8] + a[5]*b[9] + a[ 9]*b[10] + a[13]*b[11];
	c[10] = a[2]*b[8] + a[6]*b[9] + a[10]*b[10] + a[14]*b[11];
	c[11] = a[3]*b[8] + a[7]*b[9] + a[11]*b[10] + a[15]*b[11];
	c[12] = a[0]*b[12] + a[4]*b[13] + a[ 8]*b[14] + a[12]*b[15];
	c[13] = a[1]*b[12] + a[5]*b[13] + a[ 9]*b[14] + a[13]*b[15];
	c[14] = a[2]*b[12] + a[6]*b[13] + a[10]*b[14] + a[14]*b[15];
	c[15] = a[3]*b[12] + a[7]*b[13] + a[11]*b[14] + a[15]*b[15];
}

void geom_matinv2f(float  m[4]){
	float d = 1.f/(m[0]*m[3] - m[1]*m[2]);
	float t = m[0];
	m[0] = d*m[3];
	m[1] = -d*m[1];
	m[2] = -d*m[2];
	m[3] = d*t;
}
void geom_matinv2d(double m[4]){
	double d = 1./(m[0]*m[3] - m[1]*m[2]);
	double t = m[0];
	m[0] = d*m[3];
	m[1] = -d*m[1];
	m[2] = -d*m[2];
	m[3] = d*t;
}
void geom_matinv3f(float  m[9]){
	float a[9]; memcpy(a, m, 9*sizeof(float));
	float d = 1.f/(
		+a[0]*(a[4]*a[8]-a[5]*a[7])
		-a[3]*(a[1]*a[8]-a[7]*a[2])
		+a[6]*(a[1]*a[5]-a[4]*a[2])
	);
	m[0] =  (a[4]*a[8]-a[5]*a[7])*d;
	m[1] = -(a[3]*a[8]-a[6]*a[5])*d;
	m[2] =  (a[3]*a[7]-a[6]*a[4])*d;
	m[3] = -(a[1]*a[8]-a[7]*a[2])*d;
	m[4] =  (a[0]*a[8]-a[6]*a[2])*d;
	m[5] = -(a[0]*a[7]-a[1]*a[6])*d;
	m[6] =  (a[1]*a[5]-a[2]*a[4])*d;
	m[7] = -(a[0]*a[5]-a[2]*a[3])*d;
	m[8] =  (a[0]*a[4]-a[1]*a[3])*d;
}
void geom_matinv3d(double m[9]){
	double a[9]; memcpy(a, m, 9*sizeof(double));
	double d = 1./(
		+a[0]*(a[4]*a[8]-a[5]*a[7])
		-a[3]*(a[1]*a[8]-a[7]*a[2])
		+a[6]*(a[1]*a[5]-a[4]*a[2])
	);
	m[0] =  (a[4]*a[8]-a[5]*a[7])*d;
	m[1] = -(a[3]*a[8]-a[6]*a[5])*d;
	m[2] =  (a[3]*a[7]-a[6]*a[4])*d;
	m[3] = -(a[1]*a[8]-a[7]*a[2])*d;
	m[4] =  (a[0]*a[8]-a[6]*a[2])*d;
	m[5] = -(a[0]*a[7]-a[1]*a[6])*d;
	m[6] =  (a[1]*a[5]-a[2]*a[4])*d;
	m[7] = -(a[0]*a[5]-a[2]*a[3])*d;
	m[8] =  (a[0]*a[4]-a[1]*a[3])*d;
}
void geom_matinv4f(float  m[16]){
	float sum;
	unsigned int i, j, k;
	//Inversion by LU decomposition
	for(i = 1; i < 4; ++i){
		m[0+4*i] /= m[0+4*0];
	}

	for(i = 1; i < 4; ++i){
		for(j = i; j < 4; ++j){
			sum = 0.f;
			for(k = 0; k < i; ++k){
				sum += m[j+4*k] * m[k+4*i];
			}
			m[j+4*i] -= sum;
		}
		if(i == 4-1) continue;
		for(j = i+1; j < 4; ++j){
			sum = 0.f;
			for(k = 0; k < i; ++k){
				sum += m[i+4*k]*m[k+4*j];
			}
			m[i+4*j] = (m[i+4*j]-sum) / m[i+4*i];
		}
	}

	for(i = 0; i < 4; ++i){
		for(j = i; j < 4; ++j){
			sum = 1.f;
			if(i != j){
				x = 0.f;
				for(k = i; k < j; ++k){
					x -= m[j+4*k]*m[k+4*i];
				}
			}
			m[j+4*i] = x / m[j+4*j];
		}
	}
	for(i = 0; i < 4; ++i){
		for(j = i; j < 4; ++j){
			if(i == j){ continue; }
			sum = 0.f;
			for(k = i; k < j; ++k){
				sum += m[k+4*j]*( (i==k) ? 1.f : m[i+4*k] );
			}
			m[i+4*j] = -sum;
		}
	}
	for(i = 0; i < 4; ++i){
		for(j = 0; j < 4; ++j){
			sum = 0.f;
			for(k = ((i>j)?i:j); k < 4; k++ ){
				sum += ((j==k) ? 1.f : m[j+4*k])*m[k+4*i];
			}
			m[j+4*i] = sum;
		}
	}
}
void geom_matinv4d(double m[16]){
	float sum;
	unsigned int i, j, k;
	//Inversion by LU decomposition
	for(i = 1; i < 4; ++i){
		m[0+4*i] /= m[0+4*0];
	}

	for(i = 1; i < 4; ++i){
		for(j = i; j < 4; ++j){
			sum = 0.;
			for(k = 0; k < i; ++k){
				sum += m[j+4*k] * m[k+4*i];
			}
			m[j+4*i] -= sum;
		}
		if(i == 4-1) continue;
		for(j = i+1; j < 4; ++j){
			sum = 0.;
			for(k = 0; k < i; ++k){
				sum += m[i+4*k]*m[k+4*j];
			}
			m[i+4*j] = (m[i+4*j]-sum) / m[i+4*i];
		}
	}

	for(i = 0; i < 4; ++i){
		for(j = i; j < 4; ++j){
			sum = 1.;
			if(i != j){
				x = 0.;
				for(k = i; k < j; ++k){
					x -= m[j+4*k]*m[k+4*i];
				}
			}
			m[j+4*i] = x / m[j+4*j];
		}
	}
	for(i = 0; i < 4; ++i){
		for(j = i; j < 4; ++j){
			if(i == j){ continue; }
			sum = 0.;
			for(k = i; k < j; ++k){
				sum += m[k+4*j]*( (i==k) ? 1. : m[i+4*k] );
			}
			m[i+4*j] = -sum;
		}
	}
	for(i = 0; i < 4; ++i){
		for(j = 0; j < 4; ++j){
			sum = 0.;
			for(k = ((i>j)?i:j); k < 4; k++ ){
				sum += ((j==k) ? 1. : m[j+4*k])*m[k+4*i];
			}
			m[j+4*i] = sum;
		}
	}
}
