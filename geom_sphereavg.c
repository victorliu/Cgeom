#include <math.h>
#include <Cgeom/geom_la.h>

void geom_sphereavg3d(int n, const double *v, const double *w, double tol, double *avg){
	int i;
	double xVec[3] = {0,0,0};
	for(i = 0; i < n; ++i){
		xVec[0] += w[i] * v[3*i+0];
		xVec[1] += w[i] * v[3*i+1];
		xVec[2] += w[i] * v[3*i+2];
	}
	do{
		double Clocalx[3], Clocaly[3];
		// Get local basis for estimate c = xVec.:
		geom_maketriad3d(xVec, Clocalx, Clocaly);
		
		// For each v vector, compute the tangent vector from xVec 
		//		towards the v vector -- its length is the spherical length 
		//		from xVec to the v vector.  
		//	Then compute its contribution to the Hessian
		double gradient[2] = {0,0};
		double Hessian[4] = {0,0,0,0};
		for(i = 0; i < n; ++i){
			// project v[i] onto plane perpendicular to xVec, returned in vPerp
			double vPerp[3];
			{
				vPerp[0] = v[3*i+0] - xVec[0];
				vPerp[1] = v[3*i+1] - xVec[1];
				vPerp[2] = v[3*i+2] - xVec[2];
				double d = geom_dot3d(vPerp, xVec);
				vPerp[0] -= d*xVec[0];
				vPerp[1] -= d*xVec[1];
				vPerp[2] -= d*xVec[2];
			}
			double sintheta = geom_norm3d(vPerp);
			if(0. == sintheta){
				Hessian[0] += w[i];
				Hessian[3] += w[i];
			}else{
				double costheta = geom_dot3d(&v[3*i], xVec);
				double theta = atan2(sintheta, costheta); // Angle from xVec to v[i]
				double sinthetaInv = 1.0/sintheta;
				vPerp[0] *= sinthetaInv; // Normalize
				vPerp[1] *= sinthetaInv;
				vPerp[2] *= sinthetaInv;
				double cosphi = geom_dot3d(vPerp, Clocalx);
				double sinphi = geom_dot3d(vPerp, Clocaly);
				double gradlocal[2] = { cosphi, sinphi };
				gradient[0] += (w[i]*theta)*gradlocal[0]; // Added on weighted discrepancy to gradient
				gradient[1] += (w[i]*theta)*gradlocal[1];
				//Hlocal.Set( cosphi, sinphi, -sinphi, cosphi);
				//Hlocal *= LinearMapR2( w[i], 0.0, 0.0, w[i]*(theta*sinthetaInv)*costheta );
				//Hlocal *= LinearMapR2( cosphi, -sinphi, sinphi, cosphi);
				double sinphiSq = sinphi*sinphi;
				double cosphiSq = cosphi*cosphi;
				double tt = w[i]*(theta*sinthetaInv)*costheta;
				double offdiag = cosphi*sinphi*(w[i]-tt);
				Hessian[0] += w[i]*cosphiSq+tt*sinphiSq;
				Hessian[1] += offdiag;
				Hessian[2] += offdiag;
				Hessian[3] += w[i]*sinphiSq+tt*cosphiSq;
			}
		}

		double xDisplocal[2];
		geom_matinv2d(Hessian);
		geom_matvec2d(Hessian, gradient, xDisplocal);
		double xDisp[3] = {
			xDisplocal[0]*Clocalx[0] + xDisplocal[1]*Clocaly[0],
			xDisplocal[0]*Clocalx[1] + xDisplocal[1]*Clocaly[1],
			xDisplocal[0]*Clocalx[2] + xDisplocal[1]*Clocaly[2]
		};

		// cout << "    xDisp = " << xDisp << "\n";  // DEBUG

		// Step 2c: rotate xVec in direction xDisp, for new estimate.

		double xVecOld[3] = { xVec[0], xVec[1], xVec[2] };
		// Rotate unit vector xVec in the direction of xDisp: length of xDisp is rotation angle.
		//		xVec must be a unit vector.  dir must be perpindicular to xVec.
		{	
			double theta = geom_dot2d(xDisp, xDisp);
			if(0. != theta){
				theta = sqrt(theta);
				double costheta = cos(theta);
				double sintheta = sin(theta);
				double dirUnit[3] = {
					xDisp[0]/theta,
					xDisp[1]/theta,
					xDisp[2]/theta
				};
				xVec[0] = costheta*xVec[0] + sintheta*dirUnit[0];
				xVec[1] = costheta*xVec[1] + sintheta*dirUnit[1];
				xVec[2] = costheta*xVec[2] + sintheta*dirUnit[2];
			}
		}
		geom_normalize3d(xVec); // Avoid roundoff error problems

		// cout << xVec << "\n";	// DEBUG

		double diff[3] = {
			xVec[0] - xVecOld[0],
			xVec[1] - xVecOld[1],
			xVec[2] - xVecOld[2]
		};
		double errorAmt = diff[geom_imax3d(diff)];
		if(errorAmt <= tol){
			break; // return xVec as answer
		}
	}while(1);
	avg[0] = xVec[0];
	avg[1] = xVec[1];
	avg[2] = xVec[2];
}


/* returns sin(a*x)/sin(x) for -1 <= a <= 1 */
double geom_sin_ratio(double a, double x){
	double sign = 1;
	if(0 == x){ return a; }
	else if(x < 0){ x = -x; }
	
	if(0 == a){ return 0; }
	else if(a < 0){ a = -a; sign = -1; }
	if(1 == a){ return sign; }
	if(M_PI_2 == x){ return sign*sin(a*M_PI_2); }
	if(x < 1e-8){ return sign*a; }
	{
		const double x2 = x*x;
		const double a2 = a*a;
		if(x < 1e-4){ return sign*a*(1. + (1./6.)*(1.-a2)*x2); }
		else{
			return sign*sin(a*x) / sin(x);
		}
	}
}

double geom_asin_ratio(double rs, double sx){
	/* this function is odd in rs, and even in sx */
	double sign = 1;
	if(0 == sx){ return rs; }
	else if(sx < 0){ sx = -sx; }
	
	if(0 == rs){ return 0; }
	else if(rs < 0){ rs = -rs; sign = -1; }
	if(1 == rs){ return sign; }
	if(1 == sx){ return sign*asin(rs); }
	if(sx < 1e-8){ return sign*rs; }
	{
		const double s2 = sx*sx;
		const double r2 = rs*rs;
		if(sx < 1e-4){ return sign*rs*(1. - (1./6.)*(1.-r2)*s2); }
		else{
			return sign*asin(rs*sx) / asin(sx);
		}
	}
}

void geom_slerp2d(const double a[2], const double b[2], double s, double c[2]){
	if(0 == s){
		c[0] = a[0]; c[1] = a[1];
	}else if(1 == s){
		c[0] = b[0]; c[1] = b[1];
	}else{
		const double a2    = a[0]*a[0] + a[1]*a[1];
		const double dot   = a[0]*b[0] + a[1]*b[1];
		const double cross = a[0]*b[1] - a[1]*b[0];
		double theta0;
		if(fabs(dot) > fabs(cross)){ /* small or large angles */
			const double sintheta0 = cross / a2;
			theta0 = asin(sintheta0);
		}else{ /* around 90 degree angles */
			const double costheta0 = dot / a2;
			theta0 = acos(costheta0);
		}
		{
			double sr[2];
			sr[1] = geom_sin_ratio(s, theta0);
			sr[0] = geom_sin_ratio(1.-s, theta0);
			c[0] = sr[0]*a[0] + sr[1]*b[0];
			c[1] = sr[0]*a[1] + sr[1]*b[1];
		}
	}
}

double geom_unslerp2d(const double a[2], const double b[2], const double c[2]){
	const double a2    = a[0]*a[0] + a[1]*a[1];
	const double dot   = a[0]*b[0] + a[1]*b[1];
	const double cross = a[0]*b[1] - a[1]*b[0];
	const double idet = 1./cross;
	double sr[2] = {
		idet * (b[1]*c[0] - b[0]*c[1]),
		idet * (a[0]*c[1] - a[1]*c[0])
	};
	double sintheta0;
	double s[2];
	if(fabs(dot) > fabs(cross)){ /* small or large angles */
		sintheta0 = cross / a2;
	}else{ /* around 90 degree angles */
		const double costheta0 = dot / a2;
		sintheta0 = sqrt((1-costheta0)*(1+costheta0));
	}
	/* Given theta0, and sr[0] = sin(  s  *theta0)/sin(theta0),
	 *                   sr[1] = sin((1-s)*theta0)/sin(theta0),  
	 */
	s[0] = geom_asin_ratio(sr[0], sintheta0);
	s[1] = geom_asin_ratio(sr[1], sintheta0);
	{ /* fixup in case we are not an affine combination */
		double sum = s[0]+s[1];
		s[0] /= sum; s[1] /= sum;
	}
	return s[0];
}
