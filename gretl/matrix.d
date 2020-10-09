module gretl.matrix;

public import gretl.matfunctions, gretl.vector;
import gretl.base;
import std.conv, std.stdio;
import core.stdc.stdlib;
version(r) {
	import embedr.r;
}
version(standalone) {
  import std.exception;
}
version(inline) {
	private alias enforce = embedr.r.assertR;
}

private struct matrix_info {
  int t1;
  int t2;
  char **colnames;
  char **rownames;
}

struct GretlMatrix {
  int rows;
  int cols;
  double * ptr;
  matrix_info * info;

  GretlMatrix * matptr() return {
    return &this;
  }

	// The struct that wraps the underlying array should be handled by the D garbage collector
	// gretl_matrix_alloc uses malloc and leaves management of the struct itself (not just the
	// array) to the user. Therefore we copy the struct data into a D struct and destroy the
	// malloced one.
  this(int r, int c) {
    GretlMatrix * m = gretl_matrix_alloc(r, c);
    rows = m.rows;
    cols = m.cols;
    ptr = m.ptr;
    info = m.info;
    core.stdc.stdlib.free(m);
  }
  
  // Does not include col1
  this(GretlMatrix gm, int col0, int col1) {
  	rows = gm.rows;
  	cols = col1-col0;
  	ptr = &(gm.ptr[rows*col0]);
  }
  
  this(double[] v, int c=1) {
		rows = to!int(v.length)/c;
		cols = c;
		ptr = v.ptr;
	}

	// ptr and info need to be freed before destroying the struct
	// otherwise you've got a memory leak
  void free() {
    core.stdc.stdlib.free(ptr);
    core.stdc.stdlib.free(info);
  }
  
  double opIndex(int r, int c) {
    return ptr[c*rows+r];
  }

	// multidimensional slicing
	// requires the Submatrix struct
  int[2] opSlice(int dim)(int begin, int end) {
    return [begin, end];
  }
  
  Submatrix opIndex(int[2] rr, int[2] cc) {
    return Submatrix(this, rr[0], cc[0], rr[1], cc[1]);
  }

  Submatrix opIndex(int row, int[2] cc) {
    return Submatrix(this, row, cc[0], row, cc[1]);
  }
  
  Submatrix opIndex(int[2] rr, int col) {
    return Submatrix(this, rr[0], col, rr[1], col);
  }
  
  Submatrix opIndex(int row, AllElements cc) {
		return opIndex([row, row+1], [0, cols]);
	}

  Submatrix opIndex(int[2] rr, AllElements cc) {
		return opIndex(rr, [0, cols]);
	}

  Submatrix opIndex(AllElements rr, int col) {
		return opIndex([0, rows], [col, col+1]);
	}

  Submatrix opIndex(AllElements rr, int[2] cc) {
		return opIndex([0, rows], cc);
	}

  void opIndexAssign(double v, int r, int c) {
    ptr[c*rows+r] = v;
  }

  void opIndexAssign(double v, int[2] rr, int[2] cc) {
		foreach(col; cc[0]..cc[1]) {
			foreach(row; rr[0]..rr[1]) {
				this[row, col] = v;
			}
		}
  }
  
  void opIndexAssign(double v, int row, int[2] cc) {
		opIndexAssign(v, [row, row+1], cc);
  }
  
	void opIndexAssign(double v, int[2] rr, int col) {
		opIndexAssign(v, rr, [col, col+1]);
	}

	void opIndexAssign(double v, AllElements rr, int[2] cc) {
		opIndexAssign(v, [0, this.rows], cc);
	}

	void opIndexAssign(double v, int[2] rr, AllElements cc) {
		opIndexAssign(v, rr, [0, this.cols]);
	}

	void opIndexAssign(double v, AllElements rr, int col) {
		opIndexAssign(v, [0, this.rows], [col, col+1]);
	}

	void opIndexAssign(double v, int row, AllElements cc) {
		opIndexAssign(v, [row, row+1], [0, this.cols]);
	}
	
	// Assign a matrix
	void opIndexAssign(GretlMatrix m, int[2] rr, int[2] cc) {
		enforce(m.rows == rr[1]-rr[0], "Rows do not match");
		enforce(m.cols == cc[1]-cc[0], "Columns do not match");
		foreach(col; cc[0]..cc[1]) {
			foreach(row; rr[0]..rr[1]) {
				this[row, col] = m[row-rr[0], col-cc[0]];
			}
		}
	}

	void opIndexAssign(GretlMatrix m, int row, int[2] cc) {
		opIndexAssign(m, [row, row+1], cc);
	}

	void opIndexAssign(GretlMatrix m, int[2] rr, int col) {
		opIndexAssign(m, rr, [col, col+1]);
	}

	void opIndexAssign(GretlMatrix m, AllElements rr, int[2] cc) {
		opIndexAssign(m, [0, this.rows], cc);
	}

	void opIndexAssign(GretlMatrix m, int[2] rr, AllElements cc) {
		opIndexAssign(m, rr, [0, this.cols]);
	}

	void opIndexAssign(GretlMatrix m, AllElements rr, int col) {
		opIndexAssign(m, [0, this.rows], [col, col+1]);
	}

	void opIndexAssign(GretlMatrix m, int row, AllElements cc) {
		opIndexAssign(m, [row, row+1], [0, this.cols]);
	}

	// Assign a Submatrix
	void opIndexAssign(Submatrix m, int[2] rr, int[2] cc) {
		enforce(m.subRows == rr[1]-rr[0], "Rows do not match");
		enforce(m.subCols == cc[1]-cc[0], "Columns do not match");
		foreach(col; cc[0]..cc[1]) {
			foreach(row; rr[0]..rr[1]) {
				this[row, col] = m[row-rr[0], col-cc[0]];
			}
		}
	}

	void opIndexAssign(Submatrix m, int row, int[2] cc) {
		opIndexAssign(m, [row, row+1], cc);
	}

	void opIndexAssign(Submatrix m, int[2] rr, int col) {
		opIndexAssign(m, rr, [col, col+1]);
	}

	void opIndexAssign(Submatrix m, AllElements rr, int[2] cc) {
		opIndexAssign(m, [0, this.rows], cc);
	}

	void opIndexAssign(Submatrix m, int[2] rr, AllElements cc) {
		opIndexAssign(m, rr, [0, this.cols]);
	}

	void opIndexAssign(Submatrix m, AllElements rr, int col) {
		opIndexAssign(m, [0, this.rows], [col, col+1]);
	}

	void opIndexAssign(Submatrix m, int row, AllElements cc) {
		opIndexAssign(m, [row, row+1], [0, this.cols]);
	}

  // This copies, which is the expected behavior.  
  void opAssign(GretlMatrix m) {
    enforce(this.rows == m.rows, "Number of rows is different");
    enforce(this.cols == m.cols, "Number of columns is different");
    foreach(ii; 0..rows*cols) {
      ptr[ii] = m.ptr[ii];
    }
  }
  
  // Allows constructing a matrix by sending an array of rows
  void opAssign(double[][] m) {
    enforce(this.rows == m.length, "Number of rows is different");
    foreach(row, elements; m) {
			enforce(this.cols == elements.length, "Number of columns is different");
      foreach(col, val; elements) {
				this[row.to!int, col.to!int] = val;
			}
    }
  }

	// fills by column
  void opAssign(double[] v) {
    enforce(this.rows*this.cols == to!int(v.length), "double[] has different number of elements from matrix");
    foreach(ii; 0..rows*cols) {
      ptr[ii] = v[ii];
    }
  }

	version(r) {
		void opAssign(RMatrix m) {
			enforce(m.rows == this.rows, "Number of rows is different");
			enforce(m.cols == this.cols, "Number of columns is different");
			foreach(ii; 0..m.rows*m.cols) {
				ptr[ii] = m.ptr[ii];
			}
		}
	}
	
  void opAssign(DoubleMatrix m) {
    enforce(m.rows == this.rows, "Number of rows is different");
    enforce(m.cols == this.cols, "Number of columns is different");
    foreach(ii; 0..m.rows*m.cols) {
			ptr[ii] = m.data[ii];
		}
  }

  void opAssign(double a) {
    ptr[0..this.rows*this.cols] = a;
  }
  
  // Return a DoubleMatrix because returning a GretlMatrix would mean the user
  // has to manage the memory, and these functions would be very hard to use correctly.
  // How would you do x*y + z*w without some form of GC?
  DoubleMatrix opBinary(string op)(double a) {
    static if(op == "+") {
      return matrixAddition(this, a);
    }
    static if(op == "-") {
      return matrixSubtraction(this, a);
    }
    static if(op == "*") {
      return matrixMultiplication(this, a);
    }
    static if(op == "/") {
      return matrixDivision(this, a);
    }
  }

  DoubleMatrix opBinaryRight(string op)(double a) {
    static if(op == "+") {
      return matrixAddition(this, a);
    }
    static if(op == "-") {
      return matrixSubtraction(a, this);
    }
    static if(op == "*") {
      return matrixMultiplication(this, a);
    }
    static if(op == "/") {
      return matrixDivision(a, this);
    }
  }

  DoubleMatrix opBinary(string op)(GretlMatrix m) {
    static if(op == "+") {
      return matrixAddition(this, m);
    }
    static if(op == "-") {
      return matrixSubtraction(this, m);
    }
    static if(op == "*") {
      return matrixMultiplication(this, m);
    }
  }
}

struct DoubleMatrix {
  double[] data;
  int rows;
  int cols;
  private GretlMatrix temp;

  alias mat this;

  GretlMatrix mat() {
		temp.rows = rows;
		temp.cols = cols;
  	temp.ptr = data.ptr;
  	return temp;
  }
  
  GretlMatrix * matptr() return {
		temp.rows = rows;
		temp.cols = cols;
    temp.ptr = data.ptr;
    return &temp;
  }

  double * ptr() {
    return data.ptr;
  }

  this(int r, int c=1) {
    data = new double[r*c];
    rows = r;
    cols = c;
  }
  
  version(r) {
		this(RMatrix m) {
			data = new double[m.rows*m.cols];
			rows = m.rows;
			cols = m.cols;
			foreach(ii; 0..m.rows*m.cols) {
				data[ii] = m.ptr[ii];
			}
		}
	}
	
  this(double[][] m) {
		data = new double[m.length*m[0].length];
		rows = to!int(m.length);
		cols = to!int(m[0].length);
		foreach(row, vals; m) {
			foreach(col; 0..cols) {
				data[col*rows+row.to!int] = vals[col];
			}
		}
	}
	
	this(GretlMatrix * m) {
		data = new double[m.cols*m.rows];
		rows = m.rows;
		cols = m.cols;
		foreach(row; 0..rows) {
			foreach(col; 0..cols) {
				data[col*rows+row] = m.ptr[col*rows+row];
			}
		}
	}

  this(double[] v) {
    data = v;
    rows = to!int(v.length);
    cols = 1;
  }
  
  version(r) {
		RMatrix opCast(T: RMatrix)() {
			auto result = RMatrix(rows, cols);
			result.mat = this.mat;
			return result;
		}
	}
	
  double opIndex(int r, int c) {
    enforce(r < this.rows, "First index exceeds the number of rows");
    enforce(c < this.cols, "Second index exceed the number of columns");
    return data[c*this.rows + r];
  }
  
  // Support for multidimensional indexing
  int[2] opSlice(int dim)(int begin, int end) {
    return [begin, end];
  }
  
  Submatrix opIndex(int[2] rr, int[2] cc) {
    return Submatrix(this, rr[0], cc[0], rr[1], cc[1]);
  }

  Submatrix opIndex(int row, int[2] cc) {
    return opIndex([row, row+1], cc);
  }

  Submatrix opIndex(int[2] rr, int col) {
    return opIndex(rr, [col, col+1]);
  }
  
  Submatrix opIndex(int row, AllElements cc) {
		return opIndex([row, row+1], [0, cols]);
	}

  Submatrix opIndex(int[2] rr, AllElements cc) {
		return opIndex(rr, [0, cols]);
	}

  Submatrix opIndex(AllElements rr, int col) {
		return opIndex([0, rows], [col, col+1]);
	}

  Submatrix opIndex(AllElements rr, int[2] cc) {
		return opIndex([0, rows], cc);
	}

	// Assign a double to a submatrix
	void opIndexAssign(double v, int[2] rr, int[2] cc) {
		foreach(col; cc[0]..cc[1]) {
			foreach(row; rr[0]..rr[1]) {
				this[row, col] = v;
			}
		}
	}

  void opIndexAssign(double v, int r, int c) {
		enforce(r < this.rows, "Row dimension out of bounds");
		enforce(c < this.cols, "Column dimension out of bounds");
    ptr[c*rows+r] = v;
  }

	void opIndexAssign(double v, int row, int[2] cc) {
		opIndexAssign(v, [row, row+1], cc);
	}

	void opIndexAssign(double v, int[2] rr, int col) {
		opIndexAssign(v, rr, [col, col+1]);
	}

	void opIndexAssign(double v, AllElements rr, int[2] cc) {
		opIndexAssign(v, [0, this.rows], cc);
	}

	void opIndexAssign(double v, int[2] rr, AllElements cc) {
		opIndexAssign(v, rr, [0, this.cols]);
	}

	void opIndexAssign(double v, AllElements rr, int col) {
		opIndexAssign(v, [0, this.rows], [col, col+1]);
	}

	void opIndexAssign(double v, int row, AllElements cc) {
		opIndexAssign(v, [row, row+1], [0, this.cols]);
	}
	
	// Assign a matrix
	void opIndexAssign(GretlMatrix m, int[2] rr, int[2] cc) {
		enforce(m.rows == rr[1]-rr[0], "Rows do not match");
		enforce(m.cols == cc[1]-cc[0], "Columns do not match");
		foreach(col; cc[0]..cc[1]) {
			foreach(row; rr[0]..rr[1]) {
				this[row, col] = m[row-rr[0], col-cc[0]];
			}
		}
	}

	void opIndexAssign(GretlMatrix m, int row, int[2] cc) {
		opIndexAssign(m, [row, row+1], cc);
	}

	void opIndexAssign(GretlMatrix m, int[2] rr, int col) {
		opIndexAssign(m, rr, [col, col+1]);
	}

	void opIndexAssign(GretlMatrix m, AllElements rr, int[2] cc) {
		opIndexAssign(m, [0, this.rows], cc);
	}

	void opIndexAssign(GretlMatrix m, int[2] rr, AllElements cc) {
		opIndexAssign(m, rr, [0, this.cols]);
	}

	void opIndexAssign(GretlMatrix m, AllElements rr, int col) {
		opIndexAssign(m, [0, this.rows], [col, col+1]);
	}

	void opIndexAssign(GretlMatrix m, int row, AllElements cc) {
		opIndexAssign(m, [row, row+1], [0, this.cols]);
	}

	// Assign a Submatrix
	void opIndexAssign(Submatrix m, int[2] rr, int[2] cc) {
		enforce(m.subRows == rr[1]-rr[0], "Rows do not match");
		enforce(m.subCols == cc[1]-cc[0], "Columns do not match");
		foreach(col; cc[0]..cc[1]) {
			foreach(row; rr[0]..rr[1]) {
				this[row, col] = m[row-rr[0], col-cc[0]];
			}
		}
	}

	void opIndexAssign(Submatrix m, int row, int[2] cc) {
		opIndexAssign(m, [row, row+1], cc);
	}

	void opIndexAssign(Submatrix m, int[2] rr, int col) {
		opIndexAssign(m, rr, [col, col+1]);
	}

	void opIndexAssign(Submatrix m, AllElements rr, int[2] cc) {
		opIndexAssign(m, [0, this.rows], cc);
	}

	void opIndexAssign(Submatrix m, int[2] rr, AllElements cc) {
		opIndexAssign(m, rr, [0, this.cols]);
	}

	void opIndexAssign(Submatrix m, AllElements rr, int col) {
		opIndexAssign(m, [0, this.rows], [col, col+1]);
	}

	void opIndexAssign(Submatrix m, int row, AllElements cc) {
		opIndexAssign(m, [row, row+1], [0, this.cols]);
	}

  // This copies, which is the expected behavior.  
  void opAssign(GretlMatrix m) {
    enforce(this.rows == m.rows, "Number of rows is different");
    enforce(this.cols == m.cols, "Number of columns is different");
    foreach(ii; 0..rows*cols) {
      ptr[ii] = m.ptr[ii];
    }
  }
  
  // Allows constructing a matrix by sending an array of rows
  void opAssign(double[][] m) {
    enforce(this.rows == m.length, "Number of rows is different");
    foreach(row, elements; m) {
			enforce(this.cols == elements.length, "Number of columns is different");
      foreach(col, val; elements) {
				this[row.to!int, col.to!int] = val;
			}
    }
  }

	// fills by column
  void opAssign(double[] v) {
    enforce(this.rows*this.cols == to!int(v.length), "double[] has different number of elements from matrix");
    foreach(ii; 0..rows*cols) {
      ptr[ii] = v[ii];
    }
  }

	version(r) {
		void opAssign(RMatrix m) {
			enforce(m.rows == this.rows, "Number of rows is different");
			enforce(m.cols == this.cols, "Number of columns is different");
			foreach(ii; 0..m.rows*m.cols) {
				data[ii] = m.ptr[ii];
			}
		}
	}
	
  void opAssign(double a) {
    data[0..this.rows*this.cols] = a;
  }
  
  void opAssign(DoubleMatrix m) {
		data = m.data;
    rows = m.rows;
    cols = m.cols;
  }
}

// Row, Col, Submatrix do not do any reference counting.
// Up to you to make sure the reference doesn't outlive the underlying matrix.
// They are designed to be short-lived, for convenience, not for actual data storage.
struct Submatrix {
  // Original matrix
  double * ptr;
  int rows;

  // The submatrix
  int rowOffset;
  int colOffset;
  int subRows;
  int subCols;

  DoubleMatrix dup() {
    auto result = DoubleMatrix(subRows, subCols);
    foreach(col; 0..subCols) {
      foreach(row; 0..subRows) {
        result[row, col] = this[row, col];
      }
    }
    return result;
  }

  alias dup this;

  this(GretlMatrix m, int r0, int c0, int r1, int c1) {
    ptr = m.ptr;
    rows = m.rows;
    subRows = r1-r0;
    subCols = c1-c0;
    rowOffset = r0;
    colOffset = c0;
  }

  this(GretlMatrix m) {
    ptr = m.ptr;
    rows = m.rows;
    subRows = m.rows;
    subCols = m.cols;
    rowOffset = 0;
    colOffset = 0;
  }
  
  /* It's handy to be able to convert a submatrix that has only one row
   * or one column into an array.
   */
  double[] array() {
		enforce( (subCols == 1) | (subRows == 1), "Cannot convert a submatrix with multiple rows and columns into an array");
		double[] result;
		if (subCols == 1) {
			foreach(row; 0..subRows) {
				result ~= this[row,0];
			}
		} else {
			foreach(col; 0..subCols) {
				result ~= this[0,col];
			}
		}
		return result;
	}

  double opIndex(int r, int c) {
    enforce(r < subRows, "First index on Submatrix has to be less than the number of rows");
    enforce(c < subCols, "Second index on Submatrix has to be less than the number of columns");
    int newr = r+rowOffset;
    int newc = c+colOffset;
    return ptr[newc*rows + newr];
  }

  void opIndexAssign(double v, int r, int c) {
    enforce(r < subRows, "First index on Submatrix has to be less than the number of rows");
    enforce(c < subCols, "Second index on Submatrix has to be less than the number of columns");
    int newr = r+rowOffset;
    int newc = c+colOffset;
    ptr[newc*rows + newr] = v;
  }

  void opAssign(double v) {
    foreach(col; 0..subCols) {
      foreach(row; 0..subRows) {
        this[row, col] = v;
      }
    }
  }

  void opAssign(Submatrix m) {
    enforce(m.subRows == this.subRows, "Number of rows does not match");
    enforce(m.subCols == this.subCols, "Number of columns does not match");
    foreach(col; 0..subCols) {
      foreach(row; 0..subRows) {
        this[row, col] = m[row, col];
      }
    }
  }

  void opAssign(GretlMatrix m) {
    enforce(m.rows == this.subRows, "Number of rows does not match");
    enforce(m.cols == this.subCols, "Number of columns does not match");
    foreach(col; 0..m.cols) {
      foreach(row; 0..m.rows) {
        this[row, col] = m[row, col];
      }
    }
  }

  // We have this function defined because there is some overhead to using alias this with a DoubleMatrix.
  // No such overhead with an RMatrix.
  void opAssign(DoubleMatrix m) {
    enforce(m.rows == this.subRows, "Number of rows does not match");
    enforce(m.cols == this.subCols, "Number of columns does not match");
    foreach(col; 0..m.cols) {
      foreach(row; 0..m.rows) {
        this[row, col] = m[row, col];
      }
    }
  }
  
  double[] opSlice(int i0, int i1) {
		enforce( (subCols == 1) | (subRows == 1), "Can only slice a submatrix with one row or one column. Other slicing of a Submatrix is not supported at this time.");
		double[] result;
		if (subCols == 1) {
			foreach(row; i0..i1) {
				result ~= this[row,0];
			}
		} else {
			foreach(col; i0..i1) {
				result ~= this[0,col];
			}
		}
		return result;
	}
  
	version(r) {
		RMatrix rmat() {
			auto result = RMatrix(subRows, subCols);
			foreach(col; 0..subCols) {
				foreach(row; 0..subRows) {
					result[row, col] = this[row, col];
				}
			}
			return result;
		}
	}   
}

struct Row {
  GretlMatrix mat;
  int row;
  double * data;
  int length;

  this(GretlMatrix m, int r) {
    mat = m;
    row = r;
    data = &m.ptr[r];
    length = m.cols;
  }
  
  double opIndex(int c) {
    enforce(c < length, "Row index out of bounds");
    return data[c*mat.rows];
  }
  
  void opIndexAssign(double val, int ii) {
    mat[row, ii] = val;
  }
  
  void opAssign(T)(T v) {
    enforce(this.length == v.length, "Attempting to copy into Row an object with the wrong number of elements");
    foreach(ii; 0..this.length) {
      mat[row, ii] = v[ii];
    }
  }

  void opAssign(double x) {
    foreach(ii; 0..this.length) { 
      mat[row, ii] = x;
    }
  }
  
  double[] opSlice(int i0, int i1) {
		double[] result;
		foreach(col; i0..i1) {
			result ~= this[col];
		}
		return result;
	}

  bool empty() { return length == 0; }
  double front() { return data[0]; }
  void popFront() {
    data = &data[mat.rows];
    length -= 1;
  }
}

struct Col {
  GretlMatrix mat;
  int col;
  double * data;
  int length;

  this(GretlMatrix m, int c) {
    mat = m;
    col = c;
    data = &m.ptr[m.rows*c];
    length = m.rows;
  }

  double opIndex(int r) {
    enforce(r < length, "Column index out of bounds");
    return data[r];
  }

  void opIndexAssign(double val, int ii) {
    mat[ii, col] = val;
  }

  void opAssign(T)(T v) {
    enforce(this.length == v.length, "Attempting to copy into Column an object with the wrong number of elements");
    foreach(ii; 0..to!int(this.length)) {
      mat[ii, col] = v[ii];
    }
  }
 
  void opAssign(double x) {
    foreach(ii; 0..this.length) { 
      mat[ii, col] = x;
    }
  }

  void opAssign(GretlMatrix m) {
    enforce(length == m.rows, "Wrong number of elements to copy into Column");
    enforce(m.cols == 1, "Cannot copy into a Column from a matrix with more than one column");
    foreach(ii; 0..this.length) {
      mat[ii, col] = m[ii, 0];
    }
  }

  double[] opSlice(int i0, int i1) {
		double[] result;
		foreach(row; i0..i1) {
			result ~= this[row];
		}
		return result;
	}

  bool empty() { return length == 0; }
  double front() { return data[0]; }
  void popFront() {
    data = &data[1];
    length -= 1;
  }
}

struct ByRow {
  GretlMatrix mat;
  int length;
  private int rowno = 0;
  
  this(GretlMatrix m) {
    mat = m;
    length = m.rows;
  }

  bool empty() { 
    return length == 0; 
  }
  
  Row front() { 
    return Row(mat, rowno); 
  }

  void popFront() {
    rowno += 1;
    length -= 1;
  }
}

struct ByColumn {
  GretlMatrix mat;
  int colno;
  int length;

  this(GretlMatrix m) {
		writeln("start of constructor");
    mat = m;
    colno = 0;
    length = m.cols;
    writeln("exiting constructor");
  }

  bool empty() { 
    return length == 0; 
  }
  
  Col front() { 
    return Col(mat, colno); 
  }
  
  void popFront() {
    colno += 1;
    length -= 1;
  }
}
