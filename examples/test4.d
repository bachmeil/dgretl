import gretlexp.dynmodel, gretl.matrix, gretl.random;
import std.stdio;

void main() {
	auto v = TS([1.5, 6.5, 7.5, 12.5, 13.5], [2018,1], 12);
	//~ gretlexp.dynmodel.print(v);
	
	randInit();
	auto x = TS(rnorm(600), [1, 1], 1);
	//~ gretlexp.dynmodel.print(x);
	ArmaFit fit = armaGretl(x, 0, 1);
	fit.print();
}
