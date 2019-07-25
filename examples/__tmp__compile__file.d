import core.runtime;
import embedr.r, gretl.base;

struct DllInfo;

extern(C) {
	void R_init_libfoo(DllInfo * info) {
		gretl_rand_init();
		Runtime.initialize();
	}
	
	void R_unload_libfoo(DllInfo * info) {
		Runtime.terminate();
	}
}

import gretl.random;
import std.stdio;

extern (C) {
	Robj genobs(Robj rn) {
		writeln(indexSample(rn.scalar!int));
		return RNil;
	}
}
