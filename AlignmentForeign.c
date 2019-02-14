#if __STDC_VERSION__ >= 199901L
#else
#include <malloc.h>
#endif

void listVarsArrayHistoriesReduce_u(double f, int n, int* ppkk, int* pskk, int z, int* prr, double* pmv)
{
    int j;
    int k;
    int a;

    for (j = 0; j<z; j++)
    {
	for (k = 1, a = prr[z*ppkk[0] + j]; k<n; k++)
	    a = pskk[k] * a + prr[z*ppkk[k] + j];
	pmv[a] += f;
    }
}

inline int toIndex(int n, int* svv, int* ivv)
{
    int k;
    int a;

    for (k = 1, a = ivv[0]; k<n; k++)
	a = svv[k] * a + ivv[k];
    return a;
}

inline int toIndexPerm(int n, int* ppp, int* svv, int* ivv)
{
    int k;
    int p;
    int a;

    for (k = 1, a = ivv[ppp[0]]; k<n; k++)
    {
	p = ppp[k];
	a = svv[p] * a + ivv[p];
    }
    return a;
}

inline int toIndexInsert(int u, int r, int q, int n, int* svv, int* ivv)
{
    int k;
    int a;

    for (k = 0, a = 0; k<u; k++)
	a = svv[k] * a + ivv[k];
    a = r*a + q;
    for (; k<n; k++)
	a = svv[k] * a + ivv[k];
    return a;
}

inline void incIndex(int n, int* svv, int* ivv)
{
    int k;
    int y;

    for (k = n - 1; k >= 0; k--)
    {
	y = ivv[k] + 1;
	if (y == svv[k])
	    ivv[k] = 0;
	else
	{
	    ivv[k] = y;
	    break;
	}
    }
}

void listListVarsArrayHistoryPairsPartitionIndependent_u(
    double z, int v, int n, int* svv, int m, int r,
    int* lyy, int* syy, int* pppp, double* aa1, double* aa2,
    double* bb1, double* bb2)
{
    int i;
    int j;
    int k;
    int a;
    double f;
#if __STDC_VERSION__ >= 199901L
    double x1 = 0;
    double x2 = 0;
    int* ppp[m];
    double pxx1[r];
    double pxx2[r];
    double* xx1[m];
    double* xx2[m];
    int ivv[n];
    int iyy[m];
#else
    double x1;
    double x2;
    int** ppp;
    double* pxx1;
    double* pxx2;
    double** xx1;
    double** xx2;
    int* ivv;
    int* iyy;
#endif

#if __STDC_VERSION__ >= 199901L
#else
    x1 = 0.0;
    x2 = 0.0;
    ppp = (int**)alloca(sizeof(int*)*m);
    pxx1 = (double*)alloca(sizeof(double)*r);
    pxx2 = (double*)alloca(sizeof(double)*r);
    xx1 = (double**)alloca(sizeof(double*)*m);
    xx2 = (double**)alloca(sizeof(double*)*m);
    ivv = (int*)alloca(sizeof(int)*n);
    iyy = (int*)alloca(sizeof(int)*m);
#endif

    for (k = 1, a = lyy[0], ppp[0] = pppp; k<m; k++)
    {
	ppp[k] = pppp + a;
	a += lyy[k];
    }

    for (k = 1, a = syy[0], xx1[0] = pxx1, xx2[0] = pxx2; k<m; k++)
    {
	xx1[k] = pxx1 + a;
	xx2[k] = pxx2 + a;
	a += syy[k];
    }

    for (i = 0; i<r; i++)
    {
	pxx1[i] = 0.0;
	pxx2[i] = 0.0;
    }

    for (i = 0; i<n; i++)
    {
	ivv[i] = 0;
    }

    for (j = 0; j<v; j++)
    {
	for (k = 0; k<m; k++)
	{
	    i = toIndexPerm(lyy[k], ppp[k], svv, ivv);
	    xx1[k][i] += aa1[j];
	    xx2[k][i] += aa2[j];
	}
	incIndex(n, svv, ivv);
    }

    if (z != 1)
    {
	f = 1 / z;
	for (k = 0; k<m; k++)
	{
	    a = syy[k];
	    for (i = 0; i<a; i++)
	    {
		xx1[k][i] *= f;
		xx2[k][i] *= f;
	    }
	}
    }

    for (i = 0; i<m; i++)
    {
	iyy[i] = 0;
    }

    if (z != 1)
    {
	for (j = 0; j<v; j++)
	{
	    for (k = 0, x1 = z, x2 = z; k<m; k++)
	    {
		a = iyy[k];
		x1 *= xx1[k][a];
		x2 *= xx2[k][a];
	    }
	    bb1[j] = x1;
	    bb2[j] = x2;
	    incIndex(m, syy, iyy);
	}
    }
    else
    {
	for (j = 0; j<v; j++)
	{
	    for (k = 1, a = iyy[0], x1 = xx1[0][a], x2 = xx2[0][a]; k<m; k++)
	    {
		a = iyy[k];
		x1 *= xx1[k][a];
		x2 *= xx2[k][a];
	    }
	    bb1[j] = x1;
	    bb2[j] = x2;
	    incIndex(m, syy, iyy);
	}
    }
}


#include <math.h>

#ifndef INFINITY
#define INFINITY 1e+308
#endif

inline double alngam(double xvalue)

/******************************************************************************/
/*
Purpose:
ALNGAM computes the logarithm of the gamma function.
Licensing:
This code is distributed under the GNU LGPL license.
Modified:
20 November 2010
Author:
Original FORTRAN77 version by Allan Macleod.
C version by John Burkardt.
Reference:
Allan Macleod,
Algorithm AS 245,
A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
Applied Statistics,
Volume 38, Number 2, 1989, pages 397-402.
Parameters:
Input, double XVALUE, the argument of the Gamma function.
Output, double ALNGAM, the logarithm of the gamma function of X.
*/
{
    double alr2pi = 0.918938533204673;
    double r1[9] = {
	-2.66685511495,
	-24.4387534237,
	-21.9698958928,
	11.1667541262,
	3.13060547623,
	0.607771387771,
	11.9400905721,
	31.4690115749,
	15.2346874070 };
    double r2[9] = {
	-78.3359299449,
	-142.046296688,
	137.519416416,
	78.6994924154,
	4.16438922228,
	47.0668766060,
	313.399215894,
	263.505074721,
	43.3400022514 };
    double r3[9] = {
	-2.12159572323E+05,
	2.30661510616E+05,
	2.74647644705E+04,
	-4.02621119975E+04,
	-2.29660729780E+03,
	-1.16328495004E+05,
	-1.46025937511E+05,
	-2.42357409629E+04,
	-5.70691009324E+02 };
    double r4[5] = {
	0.279195317918525,
	0.4917317610505968,
	0.0692910599291889,
	3.350343815022304,
	6.012459259764103 };
    double value;
    double x;
    double x1;
    double x2;
    double xlge = 510000.0;
    double xlgst = 1.0E+30;
    double y;

    x = xvalue;
    value = 0.0;
    /*
    Check the input.
    */
    if (x <= 0.0)
    {
	return INFINITY;
    }

    /*
    Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
    */
    if (x < 1.5)
    {
	if (x < 0.5)
	{
	    value = -log(x);
	    y = x + 1.0;
	    /*
	    Test whether X < machine epsilon.
	    */
	    if (y == 1.0)
	    {
		return value;
	    }
	}
	else
	{
	    value = 0.0;
	    y = x;
	    x = (x - 0.5) - 0.5;
	}

	value = value + x * ((((
	    r1[4] * y
	    + r1[3]) * y
	    + r1[2]) * y
	    + r1[1]) * y
	    + r1[0]) / ((((
		y
		+ r1[8]) * y
		+ r1[7]) * y
		+ r1[6]) * y
		+ r1[5]);

	return value;
    }
    /*
    Calculation for 1.5 <= X < 4.0.
    */
    if (x < 4.0)
    {
	y = (x - 1.0) - 1.0;

	value = y * ((((
	    r2[4] * x
	    + r2[3]) * x
	    + r2[2]) * x
	    + r2[1]) * x
	    + r2[0]) / ((((
		x
		+ r2[8]) * x
		+ r2[7]) * x
		+ r2[6]) * x
		+ r2[5]);
    }
    /*
    Calculation for 4.0 <= X < 12.0.
    */
    else if (x < 12.0)
    {
	value = ((((
	    r3[4] * x
	    + r3[3]) * x
	    + r3[2]) * x
	    + r3[1]) * x
	    + r3[0]) / ((((
		x
		+ r3[8]) * x
		+ r3[7]) * x
		+ r3[6]) * x
		+ r3[5]);
    }
    /*
    Calculation for 12.0 <= X.
    */
    else
    {
	y = log(x);
	value = x * (y - 1.0) - 0.5 * y + alr2pi;

	if (x <= xlge)
	{
	    x1 = 1.0 / x;
	    x2 = x1 * x1;

	    value = value + x1 * ((
		r4[2] *
		x2 + r4[1]) *
		x2 + r4[0]) / ((
		    x2 + r4[4]) *
		    x2 + r4[3]);
	}
    }

    return value;
}

double alngam_1(double xvalue)

/******************************************************************************/
/*
Purpose:
ALNGAM computes the logarithm of the gamma function.
Licensing:
This code is distributed under the GNU LGPL license.
Modified:
20 November 2010
Author:
Original FORTRAN77 version by Allan Macleod.
C version by John Burkardt.
Reference:
Allan Macleod,
Algorithm AS 245,
A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
Applied Statistics,
Volume 38, Number 2, 1989, pages 397-402.
Parameters:
Input, double XVALUE, the argument of the Gamma function.
Output, double ALNGAM, the logarithm of the gamma function of X.
*/
{
    static double alr2pi = 0.918938533204673;
    static double r1[9] = {
	-2.66685511495,
	-24.4387534237,
	-21.9698958928,
	11.1667541262,
	3.13060547623,
	0.607771387771,
	11.9400905721,
	31.4690115749,
	15.2346874070 };
    static double r2[9] = {
	-78.3359299449,
	-142.046296688,
	137.519416416,
	78.6994924154,
	4.16438922228,
	47.0668766060,
	313.399215894,
	263.505074721,
	43.3400022514 };
    static double r3[9] = {
	-2.12159572323E+05,
	2.30661510616E+05,
	2.74647644705E+04,
	-4.02621119975E+04,
	-2.29660729780E+03,
	-1.16328495004E+05,
	-1.46025937511E+05,
	-2.42357409629E+04,
	-5.70691009324E+02 };
    static double r4[5] = {
	0.279195317918525,
	0.4917317610505968,
	0.0692910599291889,
	3.350343815022304,
	6.012459259764103 };
    double value;
    double x;
    double x1;
    double x2;
    static double xlge = 510000.0;
    static double xlgst = 1.0E+30;
    double y;

    x = xvalue;
    value = 0.0;
    /*
    Check the input.
    */
    if (x <= 0.0)
    {
	return INFINITY;
    }

    /*
    Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
    */
    if (x < 1.5)
    {
	if (x < 0.5)
	{
	    value = -log(x);
	    y = x + 1.0;
	    /*
	    Test whether X < machine epsilon.
	    */
	    if (y == 1.0)
	    {
		return value;
	    }
	}
	else
	{
	    value = 0.0;
	    y = x;
	    x = (x - 0.5) - 0.5;
	}

	value = value + x * ((((
	    r1[4] * y
	    + r1[3]) * y
	    + r1[2]) * y
	    + r1[1]) * y
	    + r1[0]) / ((((
		y
		+ r1[8]) * y
		+ r1[7]) * y
		+ r1[6]) * y
		+ r1[5]);

	return value;
    }
    /*
    Calculation for 1.5 <= X < 4.0.
    */
    if (x < 4.0)
    {
	y = (x - 1.0) - 1.0;

	value = y * ((((
	    r2[4] * x
	    + r2[3]) * x
	    + r2[2]) * x
	    + r2[1]) * x
	    + r2[0]) / ((((
		x
		+ r2[8]) * x
		+ r2[7]) * x
		+ r2[6]) * x
		+ r2[5]);
    }
    /*
    Calculation for 4.0 <= X < 12.0.
    */
    else if (x < 12.0)
    {
	value = ((((
	    r3[4] * x
	    + r3[3]) * x
	    + r3[2]) * x
	    + r3[1]) * x
	    + r3[0]) / ((((
		x
		+ r3[8]) * x
		+ r3[7]) * x
		+ r3[6]) * x
		+ r3[5]);
    }
    /*
    Calculation for 12.0 <= X.
    */
    else
    {
	y = log(x);
	value = x * (y - 1.0) - 0.5 * y + alr2pi;

	if (x <= xlge)
	{
	    x1 = 1.0 / x;
	    x2 = x1 * x1;

	    value = value + x1 * ((
		r4[2] *
		x2 + r4[1]) *
		x2 + r4[0]) / ((
		    x2 + r4[4]) *
		    x2 + r4[3]);
	}
    }

    return value;
}

/******************************************************************************/

double lngamma(double z)

/******************************************************************************/
/*
Purpose:
LNGAMMA computes Log(Gamma(X)) using a Lanczos approximation.
Discussion:
This algorithm is not part of the Applied Statistics algorithms.
It is slower but gives 14 or more significant decimal digits
accuracy, except around X = 1 and X = 2.   The Lanczos series from
which this algorithm is derived is interesting in that it is a
convergent series approximation for the gamma function, whereas
the familiar series due to De Moivre (and usually wrongly called
the Stirling approximation) is only an asymptotic approximation, as
is the true and preferable approximation due to Stirling.
Licensing:
This code is distributed under the GNU LGPL license.
Modified:
20 November 2010
Author:
Original FORTRAN77 version by Alan Miller.
C version by John Burkardt.
Reference:
Cornelius Lanczos,
A precision approximation of the gamma function,
SIAM Journal on Numerical Analysis, B,
Volume 1, 1964, pages 86-96.
Parameters:
Input, double Z, the argument of the Gamma function.
Output, double LNGAMMA, the logarithm of the gamma function of Z.
*/
{
    static double a[9] = {
	0.9999999999995183,
	676.5203681218835,
	-1259.139216722289,
	771.3234287757674,
	-176.6150291498386,
	12.50734324009056,
	-0.1385710331296526,
	0.9934937113930748E-05,
	0.1659470187408462E-06 };
    int j;
    static double lnsqrt2pi = 0.9189385332046727;
    double tmp;
    double value;

    if (z <= 0.0)
    {
	return INFINITY;
    }

    value = 0.0;
    tmp = z + 7.0;
    for (j = 8; 1 <= j; j--)
    {
	value = value + a[j] / tmp;
	tmp = tmp - 1.0;
    }

    value = value + a[0];
    value = log(value) + lnsqrt2pi - (z + 6.5)
	+ (z - 0.5) * log(z + 6.5);

    return value;
}


int listVarsArrayHistoriesAlignedTop_u(
    int xmax, int omax, int n, int* svv, int m, int z1, int z2,
    int* ppww, int* phh1, double* pxx1, int* phh2, double* pxx2,
    int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)
{
    int t = 0;
#if __STDC_VERSION__ >= 199901L
    double aa[xmax];
    double* xx1[n];
    double* xx2[n];
    double zf = (double)z1;
    double f = (double)z1 / (double)z2;
#else
    double* aa;
    double** xx1;
    double** xx2;
    double zf;
    double f;
#endif
    double t1;
    double t2;
    int t3;
    int tm;
    double x1;
    double x2;
    int x3;
    double a1;
    double a2;
    double b1;
    double b2;
    int ii;
    int ij;
    int pi;
    int pj;
    int qi;
    int qj;
    int si;
    int sj;
    int u;
    int i;
    int j;
    int k;
    int a;

#if __STDC_VERSION__ >= 199901L
#else
    aa = (double*)alloca(sizeof(double)*xmax);
    xx1 = (double**)alloca(sizeof(double*)*n);
    xx2 = (double**)alloca(sizeof(double*)*n);
    zf = (double)z1;
    f = (double)z1 / (double)z2;
#endif

    *s = 0;

    for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
    {
	xx1[k] = pxx1 + a;
	xx2[k] = pxx2 + a;
	a += svv[k];
    }

    for (ii = 0; ii<m - 1; ii++)
    {
	pi = ppww[ii];
	si = svv[pi];
	for (ij = ii + 1; ij<m; ij++)
	{
	    pj = ppww[ij];
	    sj = svv[pj];
	    u = si*sj;
	    if (u <= xmax)
	    {
		(*s)++;
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		qi = z1*pi;
		qj = z1*pj;
		for (j = 0; j<z1; j++)
		    aa[sj*phh1[qi + j] + phh1[qj + j]] += 1.0;
		for (a1 = 0.0, i = 0; i<u; i++)
		    a1 += alngam(aa[i]);
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		qi = z2*pi;
		qj = z2*pj;
		for (j = 0; j<z2; j++)
		    aa[sj*phh2[qi + j] + phh2[qj + j]] += f;
		for (b1 = 0.0, i = 0; i<u; i++)
		    b1 += alngam(aa[i]);
		for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
		{
		    x1 = zf*xx1[pi][i];
		    x2 = zf*xx2[pi][i];
		    for (j = 0; j<sj; j++)
		    {
			a2 += alngam(x1*xx1[pj][j] + 1.0);
			b2 += alngam(x2*xx2[pj][j] + 1.0);
		    }
		}
		if (t<omax)
		{
		    tww1[t] = pi;
		    tww2[t] = pj;
		    ts1[t] = a1 - a2 - b1 + b2;
		    ts2[t] = b2 - b1;
		    ts3[t] = -u;
		    t++;
		    if (t == omax)
		    {
			for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
			{
			    x1 = ts1[i];
			    if (t1>x1)
			    {
				t1 = x1;
				t2 = ts2[i];
				t3 = ts3[i];
				tm = i;
			    }
			    else if (t1 == x1)
			    {
				x2 = ts2[i];
				if (t2>x2)
				{
				    t2 = x2;
				    t3 = ts3[i];
				    tm = i;
				}
				else if (t2 == x2)
				{
				    x3 = ts3[i];
				    if (t3>x3)
				    {
					t3 = x3;
					tm = i;
				    }
				}
			    }
			}
		    }
		}
		else
		{
		    x1 = a1 - a2 - b1 + b2;
		    x2 = b2 - b1;
		    x3 = -u;
		    if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
		    {
			tww1[tm] = pi;
			tww2[tm] = pj;
			ts1[tm] = x1;
			ts2[tm] = x2;
			ts3[tm] = x3;
			for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
			{
			    x1 = ts1[i];
			    if (t1>x1)
			    {
				t1 = x1;
				t2 = ts2[i];
				t3 = ts3[i];
				tm = i;
			    }
			    else if (t1 == x1)
			    {
				x2 = ts2[i];
				if (t2>x2)
				{
				    t2 = x2;
				    t3 = ts3[i];
				    tm = i;
				}
				else if (t2 == x2)
				{
				    x3 = ts3[i];
				    if (t3>x3)
				    {
					t3 = x3;
					tm = i;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return t;
}

/*

int listVarsArrayHistoriesAlignedTop_u_1(
    int xmax, int omax, int n, int* svv, int m, int z1, int z2,
    int* ppww, int* phh1, double* pxx1, int* phh2, double* pxx2, int* tww1, int* tww2)
{
    int t = 0;
    double aa[xmax];
    double* xx1[n];
    double* xx2[n];
    double zf = (double)z1;
    double f = (double)z1 / (double)z2;
    double ts1[omax];
    double ts2[omax];
    int ts3[omax];
    double t1;
    double t2;
    int t3;
    int tm;
    double x1;
    double x2;
    int x3;
    double a1;
    double a2;
    double b1;
    double b2;
    int ii;
    int ij;
    int pi;
    int pj;
    int qi;
    int qj;
    int si;
    int sj;
    int u;
    int i;
    int j;
    int k;
    int a;

    for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
    {
	xx1[k] = pxx1 + a;
	xx2[k] = pxx2 + a;
	a += svv[k];
    }

    for (ii = 0; ii<m - 1; ii++)
    {
	pi = ppww[ii];
	si = svv[pi];
	for (ij = ii + 1; ij<m; ij++)
	{
	    pj = ppww[ij];
	    sj = svv[pj];
	    u = si*sj;
	    if (u <= xmax)
	    {
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		qi = z1*pi;
		qj = z1*pj;
		for (j = 0; j<z1; j++)
		    aa[sj*phh1[qi + j] + phh1[qj + j]] += 1.0;
		for (a1 = 0.0, i = 0; i<u; i++)
		    a1 += alngam(aa[i]);
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		qi = z2*pi;
		qj = z2*pj;
		for (j = 0; j<z2; j++)
		    aa[sj*phh2[qi + j] + phh2[qj + j]] += f;
		for (b1 = 0.0, i = 0; i<u; i++)
		    b1 += alngam(aa[i]);
		for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
		{
		    x1 = zf*xx1[pi][i];
		    x2 = zf*xx2[pi][i];
		    for (j = 0; j<sj; j++)
		    {
			a2 += alngam(x1*xx1[pj][j] + 1.0);
			b2 += alngam(x2*xx2[pj][j] + 1.0);
		    }
		}
		if (t<omax)
		{
		    tww1[t] = pi;
		    tww2[t] = pj;
		    ts1[t] = a1 - a2 - b1 + b2;
		    ts2[t] = b2 - b1;
		    ts3[t] = -u;
		    t++;
		    if (t == omax)
		    {
			for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
			{
			    x1 = ts1[i];
			    if (t1>x1)
			    {
				t1 = x1;
				t2 = ts2[i];
				t3 = ts3[i];
				tm = i;
			    }
			    else if (t1 == x1)
			    {
				x2 = ts2[i];
				if (t2>x2)
				{
				    t2 = x2;
				    t3 = ts3[i];
				    tm = i;
				}
				else if (t2 == x2)
				{
				    x3 = ts3[i];
				    if (t3>x3)
				    {
					t3 = x3;
					tm = i;
				    }
				}
			    }
			}
		    }
		}
		else
		{
		    x1 = a1 - a2 - b1 + b2;
		    x2 = b2 - b1;
		    x3 = -u;
		    if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
		    {
			tww1[tm] = pi;
			tww2[tm] = pj;
			ts1[tm] = x1;
			ts2[tm] = x2;
			ts3[tm] = x3;
			for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
			{
			    x1 = ts1[i];
			    if (t1>x1)
			    {
				t1 = x1;
				t2 = ts2[i];
				t3 = ts3[i];
				tm = i;
			    }
			    else if (t1 == x1)
			    {
				x2 = ts2[i];
				if (t2>x2)
				{
				    t2 = x2;
				    t3 = ts3[i];
				    tm = i;
				}
				else if (t2 == x2)
				{
				    x3 = ts3[i];
				    if (t3>x3)
				    {
					t3 = x3;
					tm = i;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return t;
}
*/

inline int isdup(int e, int pi, int* pj, int qi, int* qj)
{
    int ok = 1;
    int i;
    int j;
    int pk;

    if (pi != qi)
    {
	for (i = 0, ok = 0; i<e; i++)
	    if (pi == qj[i])
	    {
		ok = 1;
		break;
	    }
    }
    for (i = 0; ok && i<e; i++)
    {
	pk = pj[i];
	if (pk != qi)
	{
	    for (j = 0, ok = 0; j<e; j++)
		if (pk == qj[j])
		{
		    ok = 1;
		    break;
		}
	}
    }
    return ok;
}

inline int hash(int n, int e, int pi, int* pj)
{
    int h;
    int i;
    int j;
    int mk;
    int pk;
    int pl;

    for (i = -1, h = 0, mk = -1, pk = n; i<e; i++, mk = pk, pk = n)
    {
	if (mk<pi)
	    pk = pi;
	for (j = 0; j<e; j++)
	{
	    pl = pj[j];
	    if (mk<pl && pl<pk)
		pk = pl;
	}
	h = n*h + pk;
    }
    return h;
}

int listVarsListTuplesArrayHistoriesAlignedTop_u(
    int dense,
    int xmax, int omax, int n, int* svv, int m, int d, int e,
    int z1, int z2,
    int* ppww, int* ppdd,
    int* phh1, double* pxx1, int* phh2, double* pxx2,
    int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)
{
    int t = 0;
    int findm = 0;
#if __STDC_VERSION__ >= 199901L
    int* pdd[d];
    double aa[xmax];
    double* xx1[n];
    double* xx2[n];
    int ts4[omax];
    double zf = (double)z1;
    double f = (double)z1 / (double)z2;
#else
    int** pdd;
    double* aa;
    double** xx1;
    double** xx2;
    int* ts4;
    double zf;
    double f;
#endif
    double t1;
    double t2;
    int t3;
    int tm;
    double x1;
    double x2;
    double y1;
    double y2;
    int x3;
    int x4;
    double c;
    double a1;
    double a2;
    double b1;
    double b2;
    int ii;
    int ij;
    int pi;
    int* pj;
    int pk;
    int qi;
#if __STDC_VERSION__ >= 199901L
    int qj[e];
#else
    int* qj;
#endif
    int qk;
    int si;
#if __STDC_VERSION__ >= 199901L
    int sj[e];
    int yj[e];
#else
    int* sj;
    int* yj;
#endif
    int sk;
    int yk;
    int u;
    int u1;
    int i;
    int j;
    int k;
    int a;
    int ok;

#if __STDC_VERSION__ >= 199901L
#else
    pdd = (int**)alloca(sizeof(int*)*d);
    aa = (double*)alloca(sizeof(double)*xmax);
    xx1 = (double**)alloca(sizeof(double*)*n);
    xx2 = (double**)alloca(sizeof(double*)*n);
    ts4 = (int*)alloca(sizeof(int)*omax);
    zf = (double)z1;
    f = (double)z1 / (double)z2;
    qj = (int*)alloca(sizeof(int)*e);
    sj = (int*)alloca(sizeof(int)*e);
    yj = (int*)alloca(sizeof(int)*e);
#endif

    *s = 0;

    for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
    {
	xx1[k] = pxx1 + a;
	xx2[k] = pxx2 + a;
	a += svv[k];
    }

    for (k = 0; k<d; k++)
	pdd[k] = ppdd + e*k;

    for (ii = 0; ii<m; ii++)
    {
	pi = ppww[ii];
	si = svv[pi];
	for (ij = 0; ij<d; ij++)
	{
	    pj = pdd[ij];
	    for (k = 0, ok = 1, u1 = 1; k<e; k++)
	    {
		pk = pj[k];
		if (pk == pi)
		{
		    ok = 0;
		    break;
		}
		sk = svv[pk];
		sj[k] = sk;
		u1 *= sk;
	    }
	    u = u1*si;
	    if (ok && u <= xmax)
	    {
		(*s)++;
		x4 = hash(n, e, pi, pj);
		for (i = 0; i<t; i++)
		    if (ts4[i] == x4 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
		    {
			ok = 0;
			break;
		    }
		if (!ok)
		    continue;
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		pk = pj[0];
		sk = sj[0];
		qi = z1*pi;
		qk = z1*pk;
		for (k = 1; k<e; k++)
		    qj[k] = z1*pj[k];
		for (j = 0; j<z1; j++)
		{
		    for (k = 1, a = sk*phh1[qi + j] + phh1[qk + j]; k<e; k++)
			a = sj[k] * a + phh1[qj[k] + j];
		    aa[a] += 1.0;
		}
		for (a1 = 0.0, i = 0; i<u; i++)
		    a1 += alngam(aa[i]);
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		qi = z2*pi;
		qk = z2*pk;
		for (k = 1; k<e; k++)
		    qj[k] = z2*pj[k];
		for (j = 0; j<z2; j++)
		{
		    for (k = 1, a = sk*phh2[qi + j] + phh2[qk + j]; k<e; k++)
			a = sj[k] * a + phh2[qj[k] + j];
		    aa[a] += f;
		}
		for (b1 = 0.0, i = 0; i<u; i++)
		    b1 += alngam(aa[i]);
		for (k = 0; k<e; k++)
		    yj[k] = 0;
		for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
		{
		    y1 = zf*xx1[pi][i];
		    y2 = zf*xx2[pi][i];
		    for (j = 0; j<u1; j++)
		    {
			x1 = y1;
			x2 = y2;
			for (k = 0; k<e; k++)
			{
			    pk = pj[k];
			    yk = yj[k];
			    x1 *= xx1[pk][yk];
			    x2 *= xx2[pk][yk];
			}
			a2 += alngam(x1 + 1.0);
			b2 += alngam(x2 + 1.0);
			incIndex(e, sj, yj);
		    }
		}
		if (t<omax)
		{
		    if (dense)
		    {
			c = pow((double)u, 1.0 / ((double)(e + 1)));
			x1 = (a1 - a2 - b1 + b2) / c;
			x2 = (b2 - b1) / c;
		    }
		    else
		    {
			x1 = a1 - a2 - b1 + b2;
			x2 = b2 - b1;
		    }
		    x3 = -u;
		    tww1[t] = ii;
		    tww2[t] = ij;
		    ts1[t] = x1;
		    ts2[t] = x2;
		    ts3[t] = x3;
		    ts4[t] = x4;
		    t++;
		    if (t == omax)
			findm = 1;
		}
		else
		{
		    if (dense)
		    {
			c = pow((double)u, 1.0 / ((double)(e + 1)));
			x1 = (a1 - a2 - b1 + b2) / c;
			x2 = (b2 - b1) / c;
		    }
		    else
		    {
			x1 = a1 - a2 - b1 + b2;
			x2 = b2 - b1;
		    }
		    x3 = -u;
		    if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
		    {
			tww1[tm] = ii;
			tww2[tm] = ij;
			ts1[tm] = x1;
			ts2[tm] = x2;
			ts3[tm] = x3;
			ts4[tm] = x4;
			findm = 1;
		    }
		}
		if (findm)
		{
		    for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
		    {
			x1 = ts1[i];
			if (t1>x1)
			{
			    t1 = x1;
			    t2 = ts2[i];
			    t3 = ts3[i];
			    tm = i;
			}
			else if (t1 == x1)
			{
			    x2 = ts2[i];
			    if (t2>x2)
			    {
				t2 = x2;
				t3 = ts3[i];
				tm = i;
			    }
			    else if (t2 == x2)
			    {
				x3 = ts3[i];
				if (t3>x3)
				{
				    t3 = x3;
				    tm = i;
				}
			    }
			}
		    }
		    findm = 0;
		}
	    }
	}
    }
    return t;
}

/*
#define ROUND 0.001

int listVarsListTuplesArrayHistoriesAlignedTop_u_1(
    int xmax, int omax, int n, int* svv, int m, int d, int e,
    int z1, int z2,
    int* ppww, int* ppdd,
    int* phh1, double* pxx1, int* phh2, double* pxx2,
    int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)
{
    int t = 0;
    int findm = 0;
    int* pdd[d];
    double aa[xmax];
    double* xx1[n];
    double* xx2[n];
    double zf = (double)z1;
    double f = (double)z1 / (double)z2;
    double t1;
    double t2;
    int t3;
    int tm;
    double x1;
    double x2;
    double y1;
    double y2;
    int x3;
    double a1;
    double a2;
    double b1;
    double b2;
    int ii;
    int ij;
    int pi;
    int* pj;
    int pk;
    int qi;
    int qj[e];
    int qk;
    int si;
    int sj[e];
    int yj[e];
    int sk;
    int yk;
    int u;
    int u1;
    int i;
    int j;
    int k;
    int a;
    int ok;

    for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
    {
	xx1[k] = pxx1 + a;
	xx2[k] = pxx2 + a;
	a += svv[k];
    }

    for (k = 0; k<d; k++)
	pdd[k] = ppdd + e*k;

    for (ii = 0; ii<m; ii++)
    {
	pi = ppww[ii];
	si = svv[pi];
	for (ij = 0; ij<d; ij++)
	{
	    pj = pdd[ij];
	    for (k = 0, ok = 1, u1 = 1; k<e; k++)
	    {
		pk = pj[k];
		if (pk == pi)
		{
		    ok = 0;
		    break;
		}
		sk = svv[pk];
		sj[k] = sk;
		u1 *= sk;
	    }
	    u = u1*si;
	    if (ok && u <= xmax)
	    {
		(*s)++;
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		pk = pj[0];
		sk = sj[0];
		qi = z1*pi;
		qk = z1*pk;
		for (k = 1; k<e; k++)
		    qj[k] = z1*pj[k];
		for (j = 0; j<z1; j++)
		{
		    for (k = 1, a = sk*phh1[qi + j] + phh1[qk + j]; k<e; k++)
			a = sj[k] * a + phh1[qj[k] + j];
		    aa[a] += 1.0;
		}
		for (a1 = 0.0, i = 0; i<u; i++)
		    a1 += alngam(aa[i]);
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		qi = z2*pi;
		qk = z2*pk;
		for (k = 1; k<e; k++)
		    qj[k] = z2*pj[k];
		for (j = 0; j<z2; j++)
		{
		    for (k = 1, a = sk*phh2[qi + j] + phh2[qk + j]; k<e; k++)
			a = sj[k] * a + phh2[qj[k] + j];
		    aa[a] += f;
		}
		for (b1 = 0.0, i = 0; i<u; i++)
		    b1 += alngam(aa[i]);
		for (k = 0; k<e; k++)
		    yj[k] = 0;
		for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
		{
		    y1 = zf*xx1[pi][i];
		    y2 = zf*xx2[pi][i];
		    for (j = 0; j<u1; j++)
		    {
			x1 = y1;
			x2 = y2;
			for (k = 0; k<e; k++)
			{
			    pk = pj[k];
			    yk = yj[k];
			    x1 *= xx1[pk][yk];
			    x2 *= xx2[pk][yk];
			}
			a2 += alngam(x1 + 1.0);
			b2 += alngam(x2 + 1.0);
			incIndex(e, sj, yj);
		    }
		}
		if (t<omax)
		{
		    x1 = a1 - a2 - b1 + b2;
		    x2 = b2 - b1;
		    x3 = -u;
		    for (i = 0, ok = 1; i<t; i++)
			if (ts1[i]<x1 + ROUND && ts1[i]>x1 - ROUND &&
			    ts2[i]<x2 + ROUND && ts2[i]>x2 - ROUND &&
			    ts3[i] == x3 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
			{
			    ok = 0;
			    break;
			}
		    if (!ok)
			continue;
		    tww1[t] = ii;
		    tww2[t] = ij;
		    ts1[t] = x1;
		    ts2[t] = x2;
		    ts3[t] = x3;
		    t++;
		    if (t == omax)
			findm = 1;
		}
		else
		{
		    x1 = a1 - a2 - b1 + b2;
		    x2 = b2 - b1;
		    x3 = -u;
		    if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
		    {
			for (i = 0, ok = 1; i<t; i++)
			    if (ts1[i]<x1 + ROUND && ts1[i]>x1 - ROUND &&
				ts2[i]<x2 + ROUND && ts2[i]>x2 - ROUND &&
				ts3[i] == x3 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
			    {
				ok = 0;
				break;
			    }
			if (!ok)
			    continue;
			tww1[tm] = ii;
			tww2[tm] = ij;
			ts1[tm] = x1;
			ts2[tm] = x2;
			ts3[tm] = x3;
			findm = 1;
		    }
		}
		if (findm)
		{
		    for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
		    {
			x1 = ts1[i];
			if (t1>x1)
			{
			    t1 = x1;
			    t2 = ts2[i];
			    t3 = ts3[i];
			    tm = i;
			}
			else if (t1 == x1)
			{
			    x2 = ts2[i];
			    if (t2>x2)
			    {
				t2 = x2;
				t3 = ts3[i];
				tm = i;
			    }
			    else if (t2 == x2)
			    {
				x3 = ts3[i];
				if (t3>x3)
				{
				    t3 = x3;
				    tm = i;
				}
			    }
			}
		    }
		    findm = 0;
		}
	    }
	}
    }
    return t;
}

int listVarsListTuplesArrayHistoriesAlignedTop_u_2(
    int dense,
    int xmax, int omax, int n, int* svv, int m, int d, int e,
    int z1, int z2,
    int* ppww, int* ppdd,
    int* phh1, double* pxx1, int* phh2, double* pxx2,
    int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)
{
    int t = 0;
    int findm = 0;
    int* pdd[d];
    double aa[xmax];
    double* xx1[n];
    double* xx2[n];
    double zf = (double)z1;
    double f = (double)z1 / (double)z2;
    double t1;
    double t2;
    int t3;
    int tm;
    double x1;
    double x2;
    double y1;
    double y2;
    int x3;
    double c;
    double a1;
    double a2;
    double b1;
    double b2;
    int ii;
    int ij;
    int pi;
    int* pj;
    int pk;
    int qi;
    int qj[e];
    int qk;
    int si;
    int sj[e];
    int yj[e];
    int sk;
    int yk;
    int u;
    int u1;
    int i;
    int j;
    int k;
    int a;
    int ok;

    for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
    {
	xx1[k] = pxx1 + a;
	xx2[k] = pxx2 + a;
	a += svv[k];
    }

    for (k = 0; k<d; k++)
	pdd[k] = ppdd + e*k;

    for (ii = 0; ii<m; ii++)
    {
	pi = ppww[ii];
	si = svv[pi];
	for (ij = 0; ij<d; ij++)
	{
	    pj = pdd[ij];
	    for (k = 0, ok = 1, u1 = 1; k<e; k++)
	    {
		pk = pj[k];
		if (pk == pi)
		{
		    ok = 0;
		    break;
		}
		sk = svv[pk];
		sj[k] = sk;
		u1 *= sk;
	    }
	    u = u1*si;
	    if (ok && u <= xmax)
	    {
		(*s)++;
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		pk = pj[0];
		sk = sj[0];
		qi = z1*pi;
		qk = z1*pk;
		for (k = 1; k<e; k++)
		    qj[k] = z1*pj[k];
		for (j = 0; j<z1; j++)
		{
		    for (k = 1, a = sk*phh1[qi + j] + phh1[qk + j]; k<e; k++)
			a = sj[k] * a + phh1[qj[k] + j];
		    aa[a] += 1.0;
		}
		for (a1 = 0.0, i = 0; i<u; i++)
		    a1 += alngam(aa[i]);
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		qi = z2*pi;
		qk = z2*pk;
		for (k = 1; k<e; k++)
		    qj[k] = z2*pj[k];
		for (j = 0; j<z2; j++)
		{
		    for (k = 1, a = sk*phh2[qi + j] + phh2[qk + j]; k<e; k++)
			a = sj[k] * a + phh2[qj[k] + j];
		    aa[a] += f;
		}
		for (b1 = 0.0, i = 0; i<u; i++)
		    b1 += alngam(aa[i]);
		for (k = 0; k<e; k++)
		    yj[k] = 0;
		for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
		{
		    y1 = zf*xx1[pi][i];
		    y2 = zf*xx2[pi][i];
		    for (j = 0; j<u1; j++)
		    {
			x1 = y1;
			x2 = y2;
			for (k = 0; k<e; k++)
			{
			    pk = pj[k];
			    yk = yj[k];
			    x1 *= xx1[pk][yk];
			    x2 *= xx2[pk][yk];
			}
			a2 += alngam(x1 + 1.0);
			b2 += alngam(x2 + 1.0);
			incIndex(e, sj, yj);
		    }
		}
		if (t<omax)
		{
		    if (dense)
		    {
			c = pow((double)u, 1.0 / ((double)(e + 1)));
			x1 = (a1 - a2 - b1 + b2) / c;
			x2 = (b2 - b1) / c;
		    }
		    else
		    {
			x1 = a1 - a2 - b1 + b2;
			x2 = b2 - b1;
		    }
		    x3 = -u;
		    for (i = 0, ok = 1; i<t; i++)
			if (ts1[i]<x1 + ROUND && ts1[i]>x1 - ROUND &&
			    ts2[i]<x2 + ROUND && ts2[i]>x2 - ROUND &&
			    ts3[i] == x3 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
			{
			    ok = 0;
			    break;
			}
		    if (!ok)
			continue;
		    tww1[t] = ii;
		    tww2[t] = ij;
		    ts1[t] = x1;
		    ts2[t] = x2;
		    ts3[t] = x3;
		    t++;
		    if (t == omax)
			findm = 1;
		}
		else
		{
		    if (dense)
		    {
			c = pow((double)u, 1.0 / ((double)(e + 1)));
			x1 = (a1 - a2 - b1 + b2) / c;
			x2 = (b2 - b1) / c;
		    }
		    else
		    {
			x1 = a1 - a2 - b1 + b2;
			x2 = b2 - b1;
		    }
		    x3 = -u;
		    if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
		    {
			for (i = 0, ok = 1; i<t; i++)
			    if (ts1[i]<x1 + ROUND && ts1[i]>x1 - ROUND &&
				ts2[i]<x2 + ROUND && ts2[i]>x2 - ROUND &&
				ts3[i] == x3 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
			    {
				ok = 0;
				break;
			    }
			if (!ok)
			    continue;
			tww1[tm] = ii;
			tww2[tm] = ij;
			ts1[tm] = x1;
			ts2[tm] = x2;
			ts3[tm] = x3;
			findm = 1;
		    }
		}
		if (findm)
		{
		    for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
		    {
			x1 = ts1[i];
			if (t1>x1)
			{
			    t1 = x1;
			    t2 = ts2[i];
			    t3 = ts3[i];
			    tm = i;
			}
			else if (t1 == x1)
			{
			    x2 = ts2[i];
			    if (t2>x2)
			    {
				t2 = x2;
				t3 = ts3[i];
				tm = i;
			    }
			    else if (t2 == x2)
			    {
				x3 = ts3[i];
				if (t3>x3)
				{
				    t3 = x3;
				    tm = i;
				}
			    }
			}
		    }
		    findm = 0;
		}
	    }
	}
    }
    return t;
}
*/


int listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u(
    int dense,
    int xmax, int omax, int n, int* svv, int m, int d, int e,
    int z1, int z2,
    int ccl, int* ppccd, int* ppccu,
    int* ppww, int* ppdd,
    int* phh1, double* pxx1, int* phh2, double* pxx2,
    int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)
{
    int t = 0;
    int findm = 0;
#if __STDC_VERSION__ >= 199901L
    int* pdd[d];
    double aa[xmax];
    double* xx1[n];
    double* xx2[n];
    int ppccx[ccl];
    int ccx;
    int ts4[omax];
    double zf = (double)z1;
    double f = (double)z1 / (double)z2;
#else
    int** pdd;
    double* aa;
    double** xx1;
    double** xx2;
    int* ppccx;
    int ccx;
    int* ts4;
    double zf;
    double f;
#endif
    double t1;
    double t2;
    int t3;
    int tm;
    double x1;
    double x2;
    double y1;
    double y2;
    int x3;
    int x4;
    double c;
    double a1;
    double a2;
    double b1;
    double b2;
    int ii;
    int ij;
    int pi;
    int* pj;
    int pk;
    int qi;
#if __STDC_VERSION__ >= 199901L
    int qj[e];
#else
    int* qj;
#endif
    int qk;
    int si;
#if __STDC_VERSION__ >= 199901L
    int sj[e];
    int yj[e];
#else
    int* sj;
    int* yj;
#endif
    int sk;
    int yk;
    int u;
    int u1;
    int i;
    int j;
    int k;
    int h;
    int a;
    int ok;

#if __STDC_VERSION__ >= 199901L
#else
    pdd = (int**)alloca(sizeof(int*)*d);
    aa = (double*)alloca(sizeof(double)*xmax);
    xx1 = (double**)alloca(sizeof(double*)*n);
    xx2 = (double**)alloca(sizeof(double*)*n);
    ppccx = (int*)alloca(sizeof(int)*ccl);
    ts4 = (int*)alloca(sizeof(int)*omax);
    zf = (double)z1;
    f = (double)z1 / (double)z2;
    qj = (int*)alloca(sizeof(int)*e);
    sj = (int*)alloca(sizeof(int)*e);
    yj = (int*)alloca(sizeof(int)*e);
#endif

    *s = 0;

    for (k = 1, a = svv[0], xx1[0] = pxx1, xx2[0] = pxx2; k<n; k++)
    {
	xx1[k] = pxx1 + a;
	xx2[k] = pxx2 + a;
	a += svv[k];
    }

    for (k = 0; k<d; k++)
	pdd[k] = ppdd + e*k;

    for (ii = 0; ii<m; ii++)
    {
	pi = ppww[ii];
	si = svv[pi];
	for (ccx = 0, h = 0; h < ccl; h++)
	    if (ppccu[h] == pi)
	    {
		ppccx[ccx] = ppccd[h];
		ccx++;
	    }
	    else if (ppccd[h] == pi)
	    {
		ppccx[ccx] = ppccu[h];
		ccx++;
	    }
	for (ij = 0; ij<d; ij++)
	{
	    pj = pdd[ij];
	    for (k = 0, ok = 1, u1 = 1; k<e; k++)
	    {
		pk = pj[k];
		if (pk == pi)
		    ok = 0;
		for (h = 0; h<ccx; h++)
		    if (pk == ppccx[h])
			ok = 0;
		if (!ok)
		    break;
		sk = svv[pk];
		sj[k] = sk;
		u1 *= sk;
	    }
	    u = u1*si;
	    if (ok && u <= xmax)
	    {
		(*s)++;
		x4 = hash(n, e, pi, pj);
		for (i = 0; i<t; i++)
		    if (ts4[i] == x4 && isdup(e, pi, pj, ppww[tww1[i]], pdd[tww2[i]]))
		    {
			ok = 0;
			break;
		    }
		if (!ok)
		    continue;
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		pk = pj[0];
		sk = sj[0];
		qi = z1*pi;
		qk = z1*pk;
		for (k = 1; k<e; k++)
		    qj[k] = z1*pj[k];
		for (j = 0; j<z1; j++)
		{
		    for (k = 1, a = sk*phh1[qi + j] + phh1[qk + j]; k<e; k++)
			a = sj[k] * a + phh1[qj[k] + j];
		    aa[a] += 1.0;
		}
		for (a1 = 0.0, i = 0; i<u; i++)
		    a1 += alngam(aa[i]);
		for (i = 0; i<u; i++)
		    aa[i] = 1.0;
		qi = z2*pi;
		qk = z2*pk;
		for (k = 1; k<e; k++)
		    qj[k] = z2*pj[k];
		for (j = 0; j<z2; j++)
		{
		    for (k = 1, a = sk*phh2[qi + j] + phh2[qk + j]; k<e; k++)
			a = sj[k] * a + phh2[qj[k] + j];
		    aa[a] += f;
		}
		for (b1 = 0.0, i = 0; i<u; i++)
		    b1 += alngam(aa[i]);
		for (k = 0; k<e; k++)
		    yj[k] = 0;
		for (a2 = 0.0, b2 = 0.0, i = 0; i<si; i++)
		{
		    y1 = zf*xx1[pi][i];
		    y2 = zf*xx2[pi][i];
		    for (j = 0; j<u1; j++)
		    {
			x1 = y1;
			x2 = y2;
			for (k = 0; k<e; k++)
			{
			    pk = pj[k];
			    yk = yj[k];
			    x1 *= xx1[pk][yk];
			    x2 *= xx2[pk][yk];
			}
			a2 += alngam(x1 + 1.0);
			b2 += alngam(x2 + 1.0);
			incIndex(e, sj, yj);
		    }
		}
		if (t<omax)
		{
		    if (dense)
		    {
			c = pow((double)u, 1.0 / ((double)(e + 1)));
			x1 = (a1 - a2 - b1 + b2) / c;
			x2 = (b2 - b1) / c;
		    }
		    else
		    {
			x1 = a1 - a2 - b1 + b2;
			x2 = b2 - b1;
		    }
		    x3 = -u;
		    tww1[t] = ii;
		    tww2[t] = ij;
		    ts1[t] = x1;
		    ts2[t] = x2;
		    ts3[t] = x3;
		    ts4[t] = x4;
		    t++;
		    if (t == omax)
			findm = 1;
		}
		else
		{
		    if (dense)
		    {
			c = pow((double)u, 1.0 / ((double)(e + 1)));
			x1 = (a1 - a2 - b1 + b2) / c;
			x2 = (b2 - b1) / c;
		    }
		    else
		    {
			x1 = a1 - a2 - b1 + b2;
			x2 = b2 - b1;
		    }
		    x3 = -u;
		    if (t1<x1 || (t1 == x1 && t2<x2) || (t1 == x1 && t2 == x2 && t3<x3))
		    {
			tww1[tm] = ii;
			tww2[tm] = ij;
			ts1[tm] = x1;
			ts2[tm] = x2;
			ts3[tm] = x3;
			ts4[tm] = x4;
			findm = 1;
		    }
		}
		if (findm)
		{
		    for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i<omax; i++)
		    {
			x1 = ts1[i];
			if (t1>x1)
			{
			    t1 = x1;
			    t2 = ts2[i];
			    t3 = ts3[i];
			    tm = i;
			}
			else if (t1 == x1)
			{
			    x2 = ts2[i];
			    if (t2>x2)
			    {
				t2 = x2;
				t3 = ts3[i];
				tm = i;
			    }
			    else if (t2 == x2)
			    {
				x3 = ts3[i];
				if (t3>x3)
				{
				    t3 = x3;
				    tm = i;
				}
			    }
			}
		    }
		    findm = 0;
		}
	    }
	}
    }
    return t;
}


int listListVarsArrayHistoryPairsSetTuplePartitionTop_u(
    int pmax, double z, int v, int n, int* svv, int q, double y1,
    int* qm, int* ql, int* qs, int* qp, double* aa1, double* aa2,
    int* tt)
{
    int t = 0;
#if __STDC_VERSION__ >= 199901L
    double bb1[v];
    double bb2[v];
    double ts1[pmax];
    double ts2[pmax];
    int ts3[pmax];
#else
    double* bb1;
    double* bb2;
    double* ts1;
    double* ts2;
    int* ts3;
#endif
    double t1;
    double t2;
    int t3;
    int tm;
    double x1;
    double x2;
    int x3;
    int p;
    int m;
    int r;
    int i;
    double a2;
    double b2;
    double c;

#if __STDC_VERSION__ >= 199901L
#else
    bb1 = (double*)alloca(sizeof(double)*v);
    bb2 = (double*)alloca(sizeof(double)*v);
    ts1 = (double*)alloca(sizeof(double)*pmax);
    ts2 = (double*)alloca(sizeof(double)*pmax);
    ts3 = (int*)alloca(sizeof(int)*pmax);
#endif

    for (p = 0; p < q; p++)
    {
	m = qm[p];
	c = pow((double)v, 1.0 / ((double)m));

	for (r = 0, i = 0; i<m; i++)
	{
	    r += (qs + n*p)[i];
	}

	listListVarsArrayHistoryPairsPartitionIndependent_u(z, v, n, svv, m, r, ql + n*p, qs + n*p, qp + n*p, aa1, aa2, bb1, bb2);

	for (a2 = 0.0, b2 = 0.0, i = 0; i<v; i++)
	{
	    a2 += alngam(bb1[i] + 1.0);
	    b2 += alngam(bb2[i] + 1.0);
	}

	if (t < pmax)
	{
	    tt[t] = p;
	    ts1[t] = (y1 - a2 + b2) / c;
	    ts2[t] = b2;
	    ts3[t] = -m;
	    t++;
	    if (t == pmax)
	    {
		for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i < pmax; i++)
		{
		    x1 = ts1[i];
		    if (t1 > x1)
		    {
			t1 = x1;
			t2 = ts2[i];
			t3 = ts3[i];
			tm = i;
		    }
		    else if (t1 == x1)
		    {
			x2 = ts2[i];
			if (t2 > x2)
			{
			    t2 = x2;
			    t3 = ts3[i];
			    tm = i;
			}
			else if (t2 == x2)
			{
			    x3 = ts3[i];
			    if (t3 > x3)
			    {
				t3 = x3;
				tm = i;
			    }
			}
		    }
		}
	    }
	}
	else
	{
	    x1 = (y1 - a2 + b2) / c;
	    x2 = b2;
	    x3 = -m;
	    if (t1 < x1 || (t1 == x1 && t2 < x2) || (t1 == x1 && t2 == x2 && t3 < x3))
	    {
		tt[tm] = p;
		ts1[tm] = x1;
		ts2[tm] = x2;
		ts3[tm] = x3;
		for (t1 = ts1[0], t2 = ts2[0], t3 = ts3[0], tm = 0, i = 1; i < pmax; i++)
		{
		    x1 = ts1[i];
		    if (t1 > x1)
		    {
			t1 = x1;
			t2 = ts2[i];
			t3 = ts3[i];
			tm = i;
		    }
		    else if (t1 == x1)
		    {
			x2 = ts2[i];
			if (t2 > x2)
			{
			    t2 = x2;
			    t3 = ts3[i];
			    tm = i;
			}
			else if (t2 == x2)
			{
			    x3 = ts3[i];
			    if (t3 > x3)
			    {
				t3 = x3;
				tm = i;
			    }
			}
		    }
		}
	    }
	}
    }
    return t;
}

int arrayHistoryPairsRollMax_u(
    int v, int n, int* svvy, int d, int nd,
    double* aay, double* aaxy, double* bby, double* bbxy,
    int* ppm)
{
    int srchd = 0;
    double fm;

#if __STDC_VERSION__ >= 199901L
    int ivv[n];
    int syy[n];
    int svv[n];
    int szz[n];
    double aa[v];
    double aax[v];
    double bb[v];
    double bbx[v];
    double aaz[v];
    double aaxz[v];
    double bbz[v];
    double bbxz[v];
    double ff[nd];
    int ppc[nd];
#else
    int* ivv;
    int* syy;
    int* svv;
    int* szz;
    double* aa;
    double* aax;
    double* bb;
    double* bbx;
    double* aaz;
    double* aaxz;
    double* bbz;
    double* bbxz;
    double* ff;
    int* ppc;
#endif
    int minv;
    int vc;
    double fc;
    int ww;
    int sw;
    int tw;
    double fw;
    int x;    
    int y;
    int r;
    int i;
    int j;
    int k;
    int q;
    int w;
    int p;
    int s;
    int is;
    int t;
    int it;
    int m;
    int u;
    double f;
    double c;

#if __STDC_VERSION__ >= 199901L
#else
    ivv = (int*)alloca(sizeof(int)*n);
    syy = (int*)alloca(sizeof(int)*n);
    svv = (int*)alloca(sizeof(int)*n);
    szz = (int*)alloca(sizeof(int)*n);
    aa = (double*)alloca(sizeof(double)*v);
    aax = (double*)alloca(sizeof(double)*v);
    bb = (double*)alloca(sizeof(double)*v);
    bbx = (double*)alloca(sizeof(double)*v);
    aaz = (double*)alloca(sizeof(double)*v);
    aaxz = (double*)alloca(sizeof(double)*v);
    bbz = (double*)alloca(sizeof(double)*v);
    bbxz = (double*)alloca(sizeof(double)*v);
    ff = (double*)alloca(sizeof(double)*nd);
    ppc = (int*)alloca(sizeof(int)*nd);
#endif

    for (j = 0; j < v; j++)
    {
	aa[j] = aay[j];
	aax[j] = aaxy[j];
	bb[j] = bby[j];
	bbx[j] = bbxy[j];
    }

    for (i = 0; i < nd; i++)
	ff[i] = 0.0;

    for (minv = 1, i = 0; i < n; i++)
    {
	minv *= 2;
        svv[i] = svvy[i];
	for (j = 0; j < d; j++)
	{
	    p = d*i + j;
	    ppm[p] = j;
	    ppc[p] = j;
	}
    }

    vc = v;

    for (i = 0; i < n; i++)
	ivv[i] = 0;
    for (fc = 0.0, i = 0; i < vc; i++)
    {
	f = alngam(aa[i] + 1.0) - alngam(aax[i] + 1.0)
	    - alngam(bb[i] + 1.0) + alngam(bbx[i] + 1.0);
	fc += f;
	for (k = 0; k < n; k++)
	    ff[d*k + ivv[k]] += f;
	incIndex(n, svv, ivv);
    }

    m = n - 1;

    for (q = 0; vc > minv; q++)
    {
	for (x = 0, w = 0; w < n; w++)
	{
	    r = svv[w];
            if (r > 2)
            {
	    	for (i = 0; i < w; i++)
		    syy[i] = svv[i];
	    	for (; i < m; i++)
		    syy[i] = svv[i+1];
	    	p = d * w;
	    	y = vc / r;
		c = 1.0/pow((double)(y*(r-1)), 1.0/((double)n));
	    	for (s = 1; s < r; s++)
		    for (t = 0; t < s; t++)
		    {
		    	for (i = 0; i < m; i++)
			    ivv[i] = 0;
		    	for (f = fc - ff[p+s] - ff[p+t], i = 0; i < y; i++)
		    	{
			    is = toIndexInsert(w, r, s, m, syy, ivv);
			    it = toIndexInsert(w, r, t, m, syy, ivv);
			    f += alngam(aa[is] + aa[it] + 1.0) - alngam(aax[is] + aax[it] + 1.0)
			        - alngam(bb[is] + bb[it] + 1.0) + alngam(bbx[is] + bbx[it] + 1.0);
			    incIndex(m, syy, ivv);
		    	}
		    	f *= c;
		    	srchd++;
		    	if ((x == 0) || (f > fw))
		    	{
                            x++;
			    ww = w;
			    sw = s;
			    tw = t;
			    fw = f;
		    	}
		    }
            }
	}
	for (i = 0; i < n; i++)
	    szz[i] = svv[i];
	r = svv[ww] - 1;
	szz[ww] = r;
	vc /= r + 1;
	vc *= r;
	for (i = 0; i < nd; i++)
            ff[i] = 0.0;
	for (i = 0; i < n; i++)
	    ivv[i] = 0;
	for (fc = 0.0, j = 0; j < vc; j++)
	{
	    u = ivv[ww];
	    if ((u < tw) || (u > tw && u < sw))
	    {
		i = toIndex(n, svv, ivv);
		aaz[j] = aa[i];
		aaxz[j] = aax[i];
		bbz[j] = bb[i];
		bbxz[j] = bbx[i];
	    }
	    else if (u == tw)
	    {
		i = toIndex(n, svv, ivv);
		aaz[j] = aa[i];
		aaxz[j] = aax[i];
		bbz[j] = bb[i];
		bbxz[j] = bbx[i];
		ivv[ww] = sw;
		i = toIndex(n, svv, ivv);
		aaz[j] += aa[i];
		aaxz[j] += aax[i];
		bbz[j] += bb[i];
		bbxz[j] += bbx[i];
		ivv[ww] = tw;
	    }
	    else
	    {
		ivv[ww]++;
		i = toIndex(n, svv, ivv);
		aaz[j] = aa[i];
		aaxz[j] = aax[i];
		bbz[j] = bb[i];
		bbxz[j] = bbx[i];
		ivv[ww]--;
	    }
	    f = alngam(aaz[j] + 1.0) - alngam(aaxz[j] + 1.0)
		- alngam(bbz[j] + 1.0) + alngam(bbxz[j] + 1.0);
	    fc += f;
	    for (k = 0; k < n; k++)
		ff[d*k + ivv[k]] += f;
	    incIndex(n, szz, ivv);
	}
	for (j = 0; j < vc; j++)
        {
	    aa[j] = aaz[j];
	    aax[j] = aaxz[j];
	    bb[j] = bbz[j];
	    bbx[j] = bbxz[j];
        }
	svv[ww] = r;
	for (j = 0; j < d; j++)
	{
	    p = d*ww + j;
	    u = ppc[p];
	    if (u == sw)
	      ppc[p] = tw;
	    else if (u > sw)
              ppc[p] = u-1;
	}
	if (q == 0 || fw > fm)
	{
	    fm = fw;
	    for (i = 0; i < nd; i++)
		ppm[i] = ppc[i];
	}
    }

    return srchd;
}


