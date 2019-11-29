import gretl.base, gretlexp.dynmodel, gretl.matrix, gretl.random;
import std.conv, std.stdio;

void main() {
	randInit();
	
	// Generate a ts variable with 125 observations
	auto x = TS(rnorm(125), [1, 1], 1);
	Dataset * d = create_new_dataset(1, 125, 0);
	scope(exit) { destroy_dataset(d); }
	writeln(*d);
	
	// Copy the data in
	foreach(ii; 0..x.length.to!int) {
		(*d)[ii,0] = x.ptr[ii];
	}
	writeln((*d)[77,0]);
	dataset_set_time_series(d, 1, 1, 1);
	GretlPrn * prn = gretl_print_new(PrintType.stdout, null);
	int[] list = [1, 0]; 
	writeln(adf_test(3, list.ptr, d, gretlopt.v, prn));
}
