// ****
// *
// * Routines from Numerical Recipes in C to perform a polynomial
// * fit to a set of data points using general linear least squares
// * fitting (Ch. 15.4) and solution by use of the normal equations
// *
// * Dr. Stephen J. Bradshaw
// *
// * Date last modified: 01/12/2016
// *
// ****


#include <stdio.h>
#include <math.h>

#include "nrutil.h"

#define TOL 1.0e-300


// ** Routine to return the values of the polynomial coefficients a[1..ma] **
void lfit(double x[], double y[], double sig[], int ndat, double a[], int ia[], int ma, double **covar, double *chisq, void (*funcs)(double, double [], int))
{
void covsrt(double **covar, int ma, int ia[], int mfit);
void gaussj(double **a, int n, double **b, int m);
int i,j,k,l,m,mfit=0;
double ym,wt,sum,sig2i,**beta,*afunc;

beta=matrix(1,ma,1,1);
afunc=vector(1,ma);
for (j=1;j<=ma;j++)
	if (ia[j]) mfit++;
if (mfit == 0) nrerror("lfit: no parameters to be fitted");
for (j=1;j<=mfit;j++) {
	for (k=1;k<=mfit;k++) covar[j][k]=0.0;
	beta[j][1]=0.0;
}
for (i=1;i<=ndat;i++) {
	(*funcs)(x[i],afunc,ma);
	ym=y[i];
	if (mfit < ma) {
		for (j=1;j<=ma;j++)
		if (!ia[j]) ym -= a[j]*afunc[j];
	}
	sig2i=1.0/SQR(sig[i]);
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			wt=afunc[l]*sig2i;
			for (j++,k=0,m=1;m<=l;m++)
				if (ia[m]) covar[j][++k] += wt*afunc[m];
			beta[j][1] += ym*wt;
		}
	}
}
for (j=2;j<=mfit;j++)
	for (k=1;k<j;k++)
		covar[k][j]=covar[j][k];
gaussj(covar,mfit,beta,1);
for (j=0,l=1;l<=ma;l++)
	if (ia[l]) a[l]=beta[++j][1];
*chisq=0.0;
for (i=1;i<=ndat;i++) {
	(*funcs)(x[i],afunc,ma);
	for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
	*chisq += SQR((y[i]-sum)/sig[i]);
}
covsrt(covar,ma,ia,mfit);
free_vector(afunc,1,ma);
free_matrix(beta,1,ma,1,1);
}


// ** Expand in storage the covariance matrix **
void covsrt(double **covar, int ma, int ia[], int mfit)
{
int i,j,k;
double swap;

for (i=mfit+1;i<=ma;i++)
	for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
k=mfit;
for (j=ma;j>=1;j--) {
	if (ia[j]) {
		for (i=1;i<=ma;i++) {swap=covar[i][k];covar[i][k]=covar[i][j];covar[i][j]=swap;}
		for (i=1;i<=ma;i++) {swap=covar[k][i];covar[k][i]=covar[j][i];covar[j][i]=swap;}
		k--;
	}
}
}


// ** Routine to solve linear equations by Gauss-Jordan elimination **
void gaussj(double **a, int n, double **b, int m)
{
int *indxc,*indxr,*ipiv;
int i,icol,irow,j,k,l,ll;
double big,dum,pivinv,temp;
double swap;

indxc=ivector(1,n);
indxr=ivector(1,n);
ipiv=ivector(1,n);
for (j=1;j<=n;j++) ipiv[j]=0;
for (i=1;i<=n;i++) {
	big=0.0;
	for (j=1;j<=n;j++)
		if (ipiv[j] != 1)
			for (k=1;k<=n;k++) {
				if (ipiv[k] == 0) {
					if (fabs(a[j][k]) >= big) {
						big=fabs(a[j][k]);
						irow=j;
						icol=k;
					}
				}
			}
	++(ipiv[icol]);
	if (irow != icol) {
		for (l=1;l<=n;l++) {swap=a[irow][l];a[irow][l]=a[icol][l];a[icol][l]=swap;}
		for (l=1;l<=m;l++) {swap=b[irow][l];b[irow][l]=b[icol][l];b[icol][l]=swap;}
	}
	indxr[i]=irow;
	indxc[i]=icol;
	if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
	pivinv=1.0/a[icol][icol];
	a[icol][icol]=1.0;
	for (l=1;l<=n;l++) a[icol][l] *= pivinv;
	for (l=1;l<=m;l++) b[icol][l] *= pivinv;
	for (ll=1;ll<=n;ll++)
		if (ll != icol) {
			dum=a[ll][icol];
			a[ll][icol]=0.0;
			for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
			for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
		}
}
for (l=n;l>=1;l--) {
	if (indxr[l] != indxc[l])
		for (k=1;k<=n;k++)
			{swap=a[k][indxr[l]];a[k][indxr[l]]=a[k][indxc[l]];a[k][indxc[l]]=swap;}
}
free_ivector(ipiv,1,n);
free_ivector(indxr,1,n);
free_ivector(indxc,1,n);
}


// ** Routine to return the values of the polynomial coefficients a[1..ma] (using Single Value Decomposition) **
void svdfit(double x[], double y[], double sig[], int ndata, double a[], int ma, double **u, double **v, double w[], double *chisq, void (*funcs)(double, double [], int))
{
void svbksb(double **u, double w[], double **v, int m, int n, double b[],double x[]);
void svdcmp(double **a, int m, int n, double w[], double **v);
int j,i;
double wmax,tmp,thresh,sum,*b,*afunc;

b=vector(1,ndata);
afunc=vector(1,ma);
for (i=1;i<=ndata;i++) {
	(*funcs)(x[i],afunc,ma);
	tmp=1.0/sig[i];
	for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
	b[i]=y[i]*tmp;
}
svdcmp(u,ndata,ma,w,v);
wmax=0.0;
for (j=1;j<=ma;j++)
	if (w[j] > wmax) wmax=w[j];
thresh=TOL*wmax;
for (j=1;j<=ma;j++)
	if (w[j] < thresh) w[j]=0.0;
svbksb(u,w,v,ndata,ma,b,a);
*chisq=0.0;
for (i=1;i<=ndata;i++) {
	(*funcs)(x[i],afunc,ma);
	for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
	*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
}
free_vector(afunc,1,ma);
free_vector(b,1,ndata);
}


// ** Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector. No input quantities are destroyed, so the routine may be called sequentially with different b’s **
void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
int jj,j,i;
double s,*tmp;

tmp=vector(1,n);
for (j=1;j<=n;j++) {
	s=0.0;
	if (w[j]) {
		for (i=1;i<=m;i++) s += u[i][j]*b[i];
		s /= w[j];
	}
	tmp[j]=s;
}
for (j=1;j<=n;j++) {
	s=0.0;
	for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
	x[j]=s;
}
free_vector(tmp,1,n);
}


// ** Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A = U ·W ·V T . The matrix U replaces a on output. The diagonal matrix of singular values W is output as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n] **
void svdcmp(double **a, int m, int n, double w[], double **v)
{
double pythag(double a, double b);
int flag,i,its,j,jj,k,l,nm;
double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

rv1=vector(1,n);
g=scale=anorm=0.0;
for (i=1;i<=n;i++) {
	l=i+1;
	rv1[i]=scale*g;
	g=s=scale=0.0;
	if (i <= m) {
		for (k=i;k<=m;k++) scale += fabs(a[k][i]);
		if (scale) {
			for (k=i;k<=m;k++) {
				a[k][i] /= scale;
				s += a[k][i]*a[k][i];
			}
			f=a[i][i];
			g = -SIGN(sqrt(s),f);
			h=f*g-s;
			a[i][i]=f-g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
				f=s/h;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (k=i;k<=m;k++) a[k][i] *= scale;
		}
	}
	w[i]=scale *g;
	g=s=scale=0.0;
	if (i <= m && i != n) {
		for (k=l;k<=n;k++) scale += fabs(a[i][k]);
		if (scale) {
			for (k=l;k<=n;k++) {
				a[i][k] /= scale;
				s += a[i][k]*a[i][k];
			}
			f=a[i][l];
			g = -SIGN(sqrt(s),f);
			h=f*g-s;
			a[i][l]=f-g;
			for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
			for (j=l;j<=m;j++) {
				for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
				for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
			}
			for (k=l;k<=n;k++) a[i][k] *= scale;
		}
	}
	anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
}
for (i=n;i>=1;i--) {
	if (i < n) {
		if (g) {
			for (j=l;j<=n;j++)
				v[j][i]=(a[i][j]/a[i][l])/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
				for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
			}
		}
		for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
	}
	v[i][i]=1.0;
	g=rv1[i];
	l=i;
}
for (i=IMIN(m,n);i>=1;i--) {
	l=i+1;
	g=w[i];
	for (j=l;j<=n;j++) a[i][j]=0.0;
	if (g) {
		g=1.0/g;
		for (j=l;j<=n;j++) {
			for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
			f=(s/a[i][i])*g;
			for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
		}
		for (j=i;j<=m;j++) a[j][i] *= g;
	} else for (j=i;j<=m;j++) a[j][i]=0.0;
	++a[i][i];
}
for (k=n;k>=1;k--) {
	for (its=1;its<=30;its++) {
		flag=1;
		for (l=k;l>=1;l--) {
			nm=l-1;
			if ((double)(fabs(rv1[l])+anorm) == anorm) {
			flag=0;
			break;
		}
		if ((double)(fabs(w[nm])+anorm) == anorm) break;
	}
	if (flag) {
		c=0.0;
		s=1.0;
		for (i=l;i<=k;i++) {
			f=s*rv1[i];
			rv1[i]=c*rv1[i];
			if ((double)(fabs(f)+anorm) == anorm) break;
			g=w[i];
			h=pythag(f,g);
			w[i]=h;
			h=1.0/h;
			c=g*h;
			s = -f*h;
			for (j=1;j<=m;j++) {
				y=a[j][nm];
				z=a[j][i];
				a[j][nm]=y*c+z*s;
				a[j][i]=z*c-y*s;
			}
		}
	}
	z=w[k];
	if (l == k) {
		if (z < 0.0) {
			w[k] = -z;
			for (j=1;j<=n;j++) v[j][k] = -v[j][k];
		}
		break;
	}
	if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
	x=w[l];
	nm=k-1;
	y=w[nm];
	g=rv1[nm];
	h=rv1[k];
	f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	g=pythag(f,1.0);
	f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	c=s=1.0;
	for (j=l;j<=nm;j++) {
		i=j+1;
		g=rv1[i];
		y=w[i];
		h=s*g;
		g=c*g;
		z=pythag(f,h);
		rv1[j]=z;
		c=f/z;
		s=h/z;
		f=x*c+g*s;
		g = g*c-x*s;
		h=y*s;
		y *= c;
		for (jj=1;jj<=n;jj++) {
			x=v[jj][j];
			z=v[jj][i];
			v[jj][j]=x*c+z*s;
			v[jj][i]=z*c-x*s;
		}
		z=pythag(f,h);
		w[j]=z;
		if (z) {
			z=1.0/z;
			c=f*z;
			s=h*z;
		}
		f=c*g+s*y;
		x=c*y-s*g;
		for (jj=1;jj<=m;jj++) {
			y=a[jj][j];
			z=a[jj][i];
			a[jj][j]=y*c+z*s;
			a[jj][i]=z*c-y*s;
		}
	}
	rv1[l]=0.0;
	rv1[k]=f;
	w[k]=x;
}
}
free_vector(rv1,1,n);
}


// ** Computes (a2 + b2)^1/2 without destructive underflow or overflow **
double pythag(double a, double b)
{
double absa,absb;

absa=fabs(a);
absb=fabs(b);
if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


// ** To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma parameters obtained by svdfit, call this routine with matrices v[1..ma][1..ma], w[1..ma] as returned from svdfit **
void svdvar(double **v, int ma, double w[], double **cvm)
{
int k,j,i;
double sum,*wti;

wti=vector(1,ma);
for (i=1;i<=ma;i++) {
	wti[i]=0.0;
	if (w[i]) wti[i]=1.0/(w[i]*w[i]);
}
for (i=1;i<=ma;i++) {
	for (j=1;j<=i;j++) {
		for (sum=0.0,k=1;k<=ma;k++) sum += v[i][k]*v[j][k]*wti[k];
		cvm[j][i]=cvm[i][j]=sum;
	}
}
free_vector(wti,1,ma);
}