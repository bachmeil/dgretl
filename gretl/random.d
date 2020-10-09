module gretl.random;

import gretl.base, gretl.matrix, gretl.matfunctions, gretl.vector;
import std.math;
version(r) {
	import embedr.r;
}
version(standalone) {
  import std.exception;
}
version(inline) {
	private alias enforce = embedr.r.assertR;
}

void setSeed(uint seed) { 
  gretl_rand_set_seed(seed);
}

uint getSeed() { 
  return gretl_rand_get_seed();
}

void randInit() {
        gretl_rand_init();
}

/* Edit out these three functions because older Gretl libraries don't have them
bool usingBM() { 
  return gretl_rand_get_box_muller() > 0;
}

void useBM() { 
  gretl_rand_set_box_muller(1);
}

void useZiggurat() { 
  gretl_rand_set_box_muller(0);
}
*/

uint genInteger(uint k) { 
  return gretl_rand_int_max(k);
}

DoubleVector runif(int n, double min=0.0, double max=1.0) {
  auto result = DoubleVector(n);
  gretl_rand_uniform_minmax(result.ptr, 0, n-1, min, max);
  return result;
}

double runif() { 
  return gretl_rand_01();
}

DoubleVector rnorm(int n, double mean=0.0, double sd=1.0) {
  auto result = DoubleVector(n);
  gretl_rand_normal_full(result.ptr, 0, n-1, mean, sd);
  return result;
}

double rnorm() { 
  return gretl_one_snormal(); 
}

DoubleVector rmvnorm(GretlMatrix mu, GretlMatrix V) {
  enforce(mu.cols == 1, "rmvnorm: mu needs to have one column");
  enforce(mu.rows == V.rows, "rmvnorm: mu and v need to have the same number of rows");
  enforce(V.rows == V.cols, "rmvnorm: v needs to be square");
  return DoubleVector(mu + chol(V)*rnorm(mu.rows));
}

DoubleVector[] rmvnorm(int n, GretlMatrix mu, GretlMatrix V) {
	DoubleVector[] result;
	result.reserve(n);
	foreach(draw; 0..n) {
		result ~= rmvnorm(mu, V);
	}
	return result;
}

int[] indexSample(int n) {
	auto result = new int[n];
	gretl_rand_int_minmax(result.ptr, n, 0, n-1);
	return result;
}

extern(C) {
	double gretl_rand_gamma_one(double shape, double scale);
	GretlMatrix * inverse_wishart_matrix(GretlMatrix * S, int v, int * err);
	int gretl_rand_F(double * a, int t1, int t2, int v1, int v2); 
}

double rf(int v1, int v2) {
	double[] result = [0.0];
	gretl_rand_F(result.ptr, 0, 1, v1, v2);
	return result[0];
}

DoubleVector rf(int n, int v1, int v2) {
	auto result = DoubleVector(n);
	gretl_rand_F(result.ptr, 0, n, v1, v2);
	return result;
}

double rgamma(double shape, double scale) {
	return gretl_rand_gamma_one(shape, scale);
}

DoubleVector rgamma(int k, double shape, double scale) {
	auto result = DoubleVector(k);
	foreach(ii; 0..k) {
		result[ii] = rgamma(shape, scale);
	}
	return result;
}

/* This avoids allocation completely.
 * Main use is when you want to reuse an array.
 * Can be used with a GretlMatrix (completely avoiding the GC),
 * RMatrix (allocated by R), or DoubleMatrix (garbage collected).
 */
void rgamma(GretlMatrix gm, double shape, double scale) {
	foreach(ii; 0..gm.rows) {
		gm[ii,0] = gretl_rand_gamma_one(shape, scale);
	}
}

/*
 * S: Scale matrix
 * v: Degrees of freedom
 */
GretlMatrix * invWishartUnsafe(GretlMatrix S, int v) {
	GretlMatrix * gm;
	int * err;
	return inverse_wishart_matrix(&S, v, err);
}

DoubleMatrix invWishart(GretlMatrix S, int v) {
	return DoubleMatrix(invWishartUnsafe(S, v));
}

/* pdf of multivariate normal with mean mu, covariance V
 * x is the vector*/
double dmvnorm(DoubleVector x, DoubleVector mu, DoubleMatrix V, bool ln=false) {
	import std.stdio: writeln;
	DoubleVector dev = x - mu;
	DoubleMatrix inv_v = inv(V);
	DoubleVector z = DoubleVector(inv_v * dev);
	double quadform = 0.0;
	foreach(val; z) {
		quadform += val^^2;
	}
	double logSqrtDetSigma = 0.0;
	foreach(ii; 0..V.rows) {
		logSqrtDetSigma += log(V[ii,ii]);
	}
	if (ln) {
		return -0.5*quadform - logSqrtDetSigma - 0.5*V.rows*log(2.0*PI);
	} else {
		return exp(-0.5*quadform - logSqrtDetSigma - 0.5*V.rows*log(2.0*PI));
	}
}
