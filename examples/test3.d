import gretl.dynmodel, gretl.matrix;

void main() {
	auto v = TS([1.5, 6.5, 7.5, 12.5, 13.5], [2018,1], 12);
	gretl.dynmodel.print(v);
}
