import gretl.random;
import std.stdio;

extern (C) {
	Robj genobs(Robj rn) {
		writeln(indexSample(rn.scalar!int));
		return RNil;
	}
}
