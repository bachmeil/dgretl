module gretl.describe;

import gretl.base, gretl.matrix;
import std.conv, std.stdio;
version(r) {
	import embedr.r;
}
version(standalone) {
  import std.exception;
}
version(inline) {
	private alias enforce = embedr.r.assertR;
}

double min(GretlMatrix m)	{
	return gretl_min(0, m.rows*m.cols-1, m.ptr);
}

double min(double[] v) {
	return gretl_min(0, v.length.to!int-1, v.ptr);
}

double max(GretlMatrix m) {
  return gretl_max(0, m.rows*m.cols-1, m.ptr);
}

double max(double[] v) {
  return gretl_max(0, v.length.to!int-1, v.ptr);
}

double sum(GretlMatrix m) {
	return gretl_sum(0, m.rows*m.cols-1, m.ptr);
}

double sum(double[] v) {
	return gretl_sum(0, v.length.to!int-1, v.ptr);
}

double mean(GretlMatrix m) {
  return gretl_mean(0, m.rows*m.cols-1, m.ptr);
}

double mean(double[] v) {
  return gretl_mean(0, v.length.to!int-1, v.ptr);
}

/*--
Send an array of data (d) and array of probabilities (p)
Overwrites the arguments -
  d is reordered
  p is the probabilities
  p is then returned
Due to overwriting arrays, copies are made first.
--*/
double[] quantiles(double[] d, double[] p) {
	double[] temp = d.dup;
	double[] result = p.dup;
	int err = gretl_array_quantiles(temp.ptr, d.length.to!int, result.ptr, p.length.to!int);
	enforce(err == 0, "Something went wrong while calculating quantiles");
	return result;
}

double[] quantiles(GretlMatrix m, double[] p) {
	DoubleMatrix temp = m.dup;
	double[] result = p.dup;
	int err = gretl_array_quantiles(temp.ptr, temp.rows*temp.cols, result.ptr, p.length.to!int);
	enforce(err == 0, "Something went wrong while calculating quantiles");
	return result;
}

/*--
Send an array of data (d) and a probability (p)
Returns a `double` with the given probability
d is reordered, so make a copy first.
--*/
double quantile(double[] d, double p) {
	double[] temp = d.dup;
	return gretl_array_quantile(temp.ptr, d.length.to!int, p);
}

double quantile(GretlMatrix m, double p) {
	DoubleMatrix temp = m.dup;
	return gretl_array_quantile(temp.ptr, temp.rows*temp.cols, p);
}

double median(GretlMatrix m) {
	return gretl_median(0, m.rows*m.cols-1, m.ptr);
}

double median(double[] x) {
	return gretl_median(0, x.length.to!int-1, x.ptr);
}

/*--
Sum of squared deviations from the mean
--*/
double sst(GretlMatrix m) {
	return gretl_sst(0, m.rows*m.cols-1, m.ptr);
}

double sst(double[] x) {
	return gretl_sst(0, x.length.to!int-1, x.ptr);
}

double var(GretlMatrix m) {
	return gretl_variance(0, m.rows*m.cols-1, m.ptr);
}

double var(double[] x) {
	return gretl_variance(0, x.length.to!int-1, x.ptr);
}

double sd(GretlMatrix m) {
	return gretl_stddev(0, m.rows*m.cols-1, m.ptr);
}

double sd(double[] x) {
	return gretl_stddev(0, x.length.to!int-1, x.ptr);
}

double cov(GretlMatrix x, GretlMatrix y) {
	enforce(x.rows == y.rows, "First argument has " ~ x.rows.to!string ~ " rows, second argument has " ~ y.rows.to!string ~ " rows.");
	enforce(x.cols == y.cols, "First argument has " ~ x.cols.to!string ~ " columnss, second argument has " ~ y.cols.to!string ~ " columns.");
	int missing;
	return gretl_covar(0, y.rows.to!int-1, x.ptr, y.ptr, &missing);
}

double cov(double[] x, double[] y) {
	enforce(y.length == x.length, "Arrays need to have the same length");
	int missing;
	return gretl_covar(0, y.length.to!int-1, x.ptr, y.ptr, &missing);
}

double cor(GretlMatrix x, GretlMatrix y) {
	enforce(x.rows == y.rows, "First argument has " ~ x.rows.to!string ~ " rows, second argument has " ~ y.rows.to!string ~ " rows.");
	enforce(x.cols == y.cols, "First argument has " ~ x.cols.to!string ~ " columnss, second argument has " ~ y.cols.to!string ~ " columns.");
	int missing;
	return gretl_corr(0, y.rows.to!int-1, x.ptr, y.ptr, &missing);
}

double cor(double[] x, double[] y) {
	enforce(y.length == x.length, "Arrays need to have the same length");
	int missing;
	return gretl_corr(0, y.length.to!int-1, x.ptr, y.ptr, &missing);
}

double skewness(GretlMatrix m) {
	return gretl_skewness(0, m.rows*m.cols-1, m.ptr);
}

double skewness(double[] x) {
	return gretl_skewness(0, x.length.to!int-1, x.ptr);
}

double kurtosis(GretlMatrix m) {
	return gretl_kurtosis(0, m.rows*m.cols-1, m.ptr);
}

double kurtosis(double[] x) {
	return gretl_kurtosis(0, x.length.to!int-1, x.ptr);
}
