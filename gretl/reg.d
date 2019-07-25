module gretl.reg;

import gretl.base, gretl.matrix, gretl.vector;
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

struct OlsFit {
  DoubleVector coef;
  DoubleMatrix vcv;
  DoubleVector resids;
  double s2;

  double sse() {
    double result = 0.0;
    foreach(e; resids.data) { 
      result += e*e; 
    }
    return result;
  }

  void print() {
    coef.mat.print("Coefficients:");
    vcv.print("Coefficient covariance matrix:");
  }
}

OlsFit lm(GretlMatrix y, GretlMatrix x) {
  enforce(y.rows == x.rows, "y and x have different number of rows");
  enforce(y.cols == 1, "y can only have one column");
  auto coef = DoubleVector(x.cols);
  auto vcv = DoubleMatrix(x.cols, x.cols);
  auto resids = DoubleVector(y.rows);
  double s2;
  gretl_matrix_ols(y.matptr, x.matptr, coef.matptr, vcv.matptr, resids.matptr, &s2);
  return OlsFit(coef, vcv, resids, s2);
}

OlsFit lmSubsample(GretlMatrix y, GretlMatrix x, int r1, int r2) {
  return lm(y[r1..r2+1, _all], x[r1..r2+1, _all]);
}

OlsFit lm(GretlMatrix y, GretlMatrix x, int col0, int col1) {
	return lm(y, GretlMatrix(x, col0, col1));
}

// This is a bit tricky if you want an intercept
// You have to put the intercept after the variable
// We will report the intercept first to be consistent
// This function is for only that case
OlsFit lm(double[] y, double[] x, int dropFirst=0) {
	enforce(x.length == 2*y.length, "y and x have different number of rows");
	auto coef = DoubleVector(2);
	auto lhs = GretlMatrix(y[dropFirst..$-1]);
	auto rhs = GretlMatrix(x[dropFirst..$-1-dropFirst],2);
	auto vcv = DoubleMatrix(2,2);
	auto resids = DoubleVector(to!int(y.length)-dropFirst);
	double s2;
	gretl_matrix_ols(lhs.matptr, rhs.matptr, coef.matptr, vcv.matptr, resids.matptr, &s2);
	double tmp = coef[1];
	coef[1] = coef[0];
	coef[0] = tmp;
	return OlsFit(coef, vcv, resids, s2);
}
