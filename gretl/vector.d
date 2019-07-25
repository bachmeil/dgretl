module gretl.vector;

import gretl.base, gretl.matrix;
import std.conv, std.stdio;
version(standalone) {
  import std.exception;
}
version(r) {
	import embedr.r;
}
// If embedding inside R, want to use R exceptions, not D exceptions.
// D exceptions cause a segfault. R exceptions exit nicely.
version(inline) {
	private alias enforce = embedr.r.assertR;
}

struct DoubleVector {
	DoubleMatrix mat;

	this(int r) {
		mat = DoubleMatrix(r,1);
	}
	
	this(DoubleMatrix m) {
		enforce(m.cols == 1, "Cannot convert a DoubleMatrix to a DoubleVector unless it has one column.");
		mat = m;
	}
	
	this(double[] v) {
		mat = DoubleMatrix(v);
	}
	
	alias mat this;
	
	double opIndex(int r) {
    enforce(r < mat.rows, "Index out of range: index on DoubleVector is too large");
    return mat[r,0];
  }

  void opIndexAssign(double v, int r) {
    enforce(r < mat.rows, "Index out of range: index on DoubleVector is too large");
    mat[r,0] = v;
  }
  
  void opIndexAssign(GretlMatrix m, int r) {
		enforce( (m.rows == 1) & (m.cols == 1), "To assign a matrix as an element of a DoubleVector, the matrix needs to have one element.");
		opIndexAssign(m[0,0], r);
	} 

  void opAssign(T)(T x) {
    enforce(x.length == mat.rows, "Cannot assign to DoubleVector from an object with the wrong length");
    foreach(ii; 0..mat.rows) {
      mat[ii,0] = x[ii];
    }
  }
	
	void opAssign(GretlMatrix gm) {
		enforce(gm.cols == 1, "Can only copy a matrix into a DoubleVector if it has one column");
		enforce(gm.rows == mat.rows, "Number of rows in matrix is different from number of rows in DoubleVector");
		foreach(ii; 0..mat.rows) {
			mat[ii,0] = gm[ii,0];
		}
	}
	
	void opAssign(double x) {
		foreach(ii; 0..mat.rows) {
			mat[ii,0] = x;
		}
	}
	
  int[2] opSlice(int dim)(int begin, int end) {
    return [begin, end];
  }

	void opIndexAssign(double v, int[2] rr) {
		mat.opIndexAssign(v, rr, 0);
	}
	
	double * ptr() {
		return mat.data.ptr;
	}
	
	int rows() {
		return mat.rows;
	}
	
	int cols() {
		return 1;
	}
	
	int length() {
		return mat.rows;
	}
	
	void print(string msg="") {
		writeln(msg);
    foreach(val; mat.data) {
      writeln(val);
    }
  }

  bool empty() {
    return mat.rows == 0;
  }

  double front() {
    return this[0];
  }

  void popFront() {
    mat.data = mat.data[1..$];
    mat.rows -= 1;
  }
}

void replace(ref DoubleVector m, DoubleVector newmat) {
	gretl.matfunctions.replace(m.mat, newmat.mat);
}

// Only use this for filling - nothing else!
struct PutArray {
	int length;
	double * ptr;
	int fillPointer = 0;

	version(r) {
		this(RVector rv) {
			length = rv.length;
			ptr = rv.ptr;
		}
	}
	
	this(int len, double * p) {
		length = len;
		ptr = p;
	}
	
	void put(double v) {
		enforce(fillPointer < length, "Attempting to put more elements into a PutArray than have been allocated. The array has only " ~ length.to!string ~ " elements allocated.");
		ptr[fillPointer] = v;
		fillPointer += 1;
	}
	
	double first() {
		enforce(fillPointer > 0, "Attempting to return the last element of an empty array");
		return ptr[0];
	}
	
	double last() {
		enforce(fillPointer > 0, "Attempting to return the last element of an empty array");
		return ptr[fillPointer-1];
	}
	
	void reset() {
		fillPointer = 0;
	}
}
