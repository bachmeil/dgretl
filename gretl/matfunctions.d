module gretl.matfunctions;

import gretl.base, gretl.matrix;
import std.algorithm.iteration, std.conv, std.range, std.stdio;
version(r) {
	import embedr.r;
}
version(standalone) {
  import std.exception;
}
version(inline) {
	private alias enforce = embedr.r.assertR;
}

static struct AllElements {};
enum AllElements _all = {};

DoubleMatrix matrixAddition(GretlMatrix m, double a) {
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows*m.cols) {
    result.data[ii] = m.ptr[ii] + a;
  }
  return result;
}

DoubleMatrix matrixSubtraction(GretlMatrix m, double a) {
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows*m.cols) {
    result.data[ii] = m.ptr[ii] - a;
  }
  return result;
}

DoubleMatrix matrixSubtraction(double a, GretlMatrix m) {
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows*m.cols) {
    result.data[ii] = a - m.ptr[ii];
  }
  return result;
}

DoubleMatrix matrixMultiplication(GretlMatrix m, double a) {
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows*m.cols) {
    result.data[ii] = a*m.ptr[ii];
  }
  return result;
}

// Access to gretl_matrix_multiply_mod
DoubleMatrix matrixMultiplyMod(GretlMatrix m1, matmod mod1, GretlMatrix m2, matmod mod2, matmod mod3=matmod.none) {
	int rows, cols;
	
	if ( (mod1 == matmod.none) & (mod2 == matmod.none) ) {
		rows = m1.rows;
		cols = m2.cols;
		enforce(m1.cols == m2.rows, "Dimensions do not allow for matrix multiplication");
	} else if ( (mod1 == matmod.transpose) & (mod2 == matmod.none) ) {
		rows = m1.cols;
		cols = m2.cols;
		enforce(m1.rows == m2.rows, "Dimensions do not allow for matrix multiplication");
	} else if ( (mod1 == matmod.none) & (mod2 == matmod.transpose) ) {
		rows = m1.rows;
		cols = m2.rows;
		enforce(m1.cols == m2.cols, "Dimensions do not allow for matrix multiplication");
	} else if ( (mod1 == matmod.transpose) & (mod2 == matmod.transpose) ) {
		rows = m1.cols;
		cols = m2.rows;
		enforce(m1.rows == m2.cols, "Dimensions do not allow for matrix multiplication");
	} else {
		enforce(false, "Invalid options passed in for mod1 and/or mod2");
	}
	
	auto result = DoubleMatrix(rows, cols);
	int err = gretl_matrix_multiply_mod(m1.matptr, mod1, m2.matptr, mod2, result.matptr, mod3);
	enforce(err == 0, "Something went wrong in the call to gretl_matrix_multiply_mod");
	return result;
}

DoubleMatrix xtx(GretlMatrix m) {
	return matrixMultiplyMod(m, matmod.transpose, m, matmod.none);
}

DoubleMatrix matrixDivision(GretlMatrix m, double a) {
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows*m.cols) {
    result.data[ii] = m.ptr[ii]/a;
  }
  return result;
}

DoubleMatrix matrixDivision(double a, GretlMatrix m) {
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows*m.cols) {
    result.data[ii] = a/m.ptr[ii];
  }
  return result;
}

DoubleMatrix matrixAddition(GretlMatrix x, GretlMatrix y) {
  enforce(x.rows == y.rows, "Differing number of rows in matrix addition");
  enforce(x.cols == y.cols, "Differing number of columns in matrix addition");
  auto result = DoubleMatrix(x.rows, x.cols);
  foreach(ii; 0..x.rows*x.cols) {
    result.data[ii] = x.ptr[ii] + y.ptr[ii];
  }
  return result;
}

DoubleMatrix matrixSubtraction(GretlMatrix x, GretlMatrix y) {
  enforce(x.rows == y.rows, "Differing number of rows in matrix subtraction");
  enforce(x.cols == y.cols, "Differing number of columns in matrix subtraction");
  auto result = DoubleMatrix(x.rows, x.cols);
  foreach(ii; 0..x.rows*x.cols) {
    result.data[ii] = x.ptr[ii] - y.ptr[ii];
  }
  return result;
}

DoubleMatrix matrixMultiplication(GretlMatrix x, GretlMatrix y) {
  enforce(x.cols == y.rows, "Dimensions do not allow for matrix multiplication");
  auto result = DoubleMatrix(x.rows, y.cols);
  gretl_matrix_multiply(x.matptr, y.matptr, result.matptr);
  return result;
}

DoubleMatrix dup(GretlMatrix m) {
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows*m.cols) {
    result.data[ii] = m.ptr[ii];
  }
  return result;
}

void print(GretlMatrix m, string msg="") {
	writeln(msg);
  foreach(row; ByRow(m)) {
    foreach(val; row) {
      write(val, " ");
    }
    writeln("");
  }
}

DoubleMatrix copyRows(GretlMatrix m, int r0, int r1) {
  auto result = DoubleMatrix(r1-r0, m.cols);
  foreach(int row; r0..r1) {
    Row(result, row-r0) = Row(m, row);
  }
  return result;
}

/* Drop r0..r1 inclusive
   Unintuitive to not include r1 */
DoubleMatrix dropRows(GretlMatrix m, int r0Raw, int r1Raw) {
  int r0, r1;
  /* If r0Raw < 0, start at 0 */
  if (r0Raw < 0) { r0 = 0; } else { r0 = r0Raw; }
  /* If r1Raw > T, end at T */
  if (r1Raw > m.rows-1) { r1 = m.rows-1; } else { r1 = r1Raw; }

  /* Special case: drop first r1 rows */
  if (r0 == 0) {
    return copyRows(m, r1+1, m.rows);
  }

  /* Special case: drop everything starting with row r0 */
  if (r1 == m.rows-1) {
    return copyRows(m, 0, r0);
  }

  /* No special case: have to combine two chunks into one matrix
     Have to copy rows 0..r0-1 and r1+1..T */
  auto result = DoubleMatrix(m.rows-(r1-r0+1), m.cols);
  foreach(int row; 0..r0) {
    Row(result, row) = Row(m, row);
  }
  foreach(int row; r1+1..m.rows) {
    Row(result, row-(r1-r0+1)) = Row(m, row);
  }
  return result;
}

version(r) {
	RMatrix dup(RMatrix rm) { 
		RMatrix result = RMatrix(Rf_protect(Rf_duplicate(rm.robj)), true);
		return result;
	}
}

void print(Row r, string msg="") {
	writeln(msg);
  foreach(val; r) {
    writeln(val);
  }
}

void print(Col c, string msg="") {
	writeln(msg);
  foreach(val; c) {
    writeln(val);
  }
}

DoubleMatrix chol(GretlMatrix m) {
  enforce(m.rows == m.cols, "You are trying to compute a Cholesky decomposition of a non-square matrix");
  auto result = dup(m);
  int err = gretl_matrix_cholesky_decomp(result.matptr);
  enforce(err == 0, "Cholesky decomposition failed");
  return result;
}

DoubleMatrix elMultiply(GretlMatrix x, GretlMatrix y) {
	enforce(x.rows == y.rows, "Number of rows do not match");
	enforce(x.cols == y.cols, "Number of columns do not match");
	
	auto result = DoubleMatrix(x.rows, x.cols);
	foreach(ii; 0..x.rows*x.cols) {
		result.data[ii] = x.ptr[ii]*y.ptr[ii];
	}
	return result;
}

DoubleMatrix inv(GretlMatrix m) {
  enforce(m.rows == m.cols, "You are trying to take the inverse of a non-square matrix");
  DoubleMatrix result = dup(m);
  int err = gretl_invert_matrix(result.matptr);
  enforce(err == 0, "Taking the inverse of a matrix failed");
  return result;
}

DoubleMatrix transpose(GretlMatrix m) {
  auto result = DoubleMatrix(m.cols, m.rows);
  int err = gretl_matrix_transpose(result.matptr, m.matptr);
  enforce(err == 0, "Taking the transpose of a matrix failed");
  return result;
}

double det(GretlMatrix m) {
  DoubleMatrix temp = dup(m);
  int err;
  double result = gretl_matrix_determinant(temp.matptr, &err);
  enforce(err == 0, "Taking the determinant of a matrix failed");
  return result;
}

double logdet(GretlMatrix m) {
  DoubleMatrix temp = dup(m);
  int err;
  double result = gretl_matrix_log_determinant(temp.matptr, &err);
  enforce(err == 0, "Taking the log determinant of a matrix failed");
  return result;
}

DoubleMatrix diag(GretlMatrix m) {
  enforce(m.rows == m.cols, "diag is not intended to be used to take the diagonal of a non-square matrix");
  auto result = DoubleMatrix(m.rows, 1);
  foreach(ii; 0..m.rows) { 
    result[ii,0] = m[ii,ii]; 
  }
  return result;
}

DoubleVector diag(string op)(GretlMatrix y, GretlMatrix x) {
  enforce(y.cols == x.rows, "Wrong dimensions for matrix multiplication");
  auto result = DoubleVector(y.rows);
  static if(op == "*") {
    //~ foreach(ii; 0..y.rows) {
      //~ result[ii] = y[ii, _all] * x[_all, ii];
    //~ }
    double ans;
    foreach(element; 0..y.rows) {
			ans = 0.0;
			foreach(ii; 0..y.cols) {
				ans += y[element, ii]*x[ii, element];
			}
			result[element] = ans;
		}
  } else {
    static assert(false, "diag currently defined only to work with * operator");
  }
  return result;
}

void setDiagonal(GretlMatrix m, GretlMatrix newdiag) {
  enforce(newdiag.cols == 1, "Cannot set a diagonal to a matrix with more than one column");
  enforce(m.rows == m.cols, "Cannot set the diagonal of a non-square matrix");
  enforce(m.rows == newdiag.rows, "Wrong number of elements in the new diagonal");
  foreach(ii; 0..newdiag.rows) {
    m[ii,ii] = newdiag[ii,0]; 
  }
}

void setDiagonal(T)(GretlMatrix m, T v) {
  enforce(m.rows == m.cols, "Cannot set the diagonal of a non-square matrix");
  enforce(m.rows == v.length, "Wrong number of elements in the new diagonal");
  foreach(ii; 0..m.rows) { 
    m[ii,ii] = v[ii]; 
  }
}

void setDiagonal(GretlMatrix m, double v) {
  enforce(m.rows == m.cols, "Attempting to set the diagonal of a non-square matrix");
  foreach(ii; 0..m.rows) { 
    m[ii,ii] = v; 
  }
}

double trace(GretlMatrix m) { 
  enforce(m.rows == m.cols, "Cannot take the trace of a non-square matrix");
  return gretl_matrix_trace(m.matptr); 
}

DoubleMatrix raise(GretlMatrix m, double k) {
  DoubleMatrix result = dup(m);
  gretl_matrix_raise(result.matptr, k);
  return result;
}

DoubleMatrix solve(GretlMatrix x, GretlMatrix y) {
  auto temp = dup(x);
  auto result = dup(y);
  int err = gretl_LU_solve(temp.matptr, result.matptr);
  enforce(err == 0, "Call to solve failed");
  return result;
}

/*
# eye

Return a `k x k` identity matrix.
*/

DoubleMatrix eye(int k) {
  auto result = DoubleMatrix(k, k);
  result = 0.0;
  foreach(ii; 0..k) { 
    result[ii,ii] = 1.0;
  }
  return result;
}

/*
# kron

Returns the kronecker product $x \otimes y$.
*/

DoubleMatrix kron(GretlMatrix x, GretlMatrix y) {
  auto result = DoubleMatrix(x.rows*y.rows, x.cols*y.cols);
  int err = gretl_matrix_kronecker_product(x.matptr, y.matptr, result.matptr);
  enforce(err == 0, "Kronecker product failed");
  return result;
}

/*
# `middleColumns`

Return columns `c0..c1` of `gm` with no copying.

**Warning:** Similar to a D array slice, this struct holds a reference 
to the data. If the original data matrix is freed, anything can happen
when you work with this matrix. The motivation for holding only a 
reference to the data is speed.
*/

GretlMatrix middleColumns(GretlMatrix gm, int c0, int c1) {
	return GretlMatrix(gm, c0, c1);
}

/*
# `firstColumns`

Return first `k` columns of `gm` with no copying.

**Warning:** See the warning on `middleColumns`.
*/

GretlMatrix firstColumns(GretlMatrix gm, int k) {
	return GretlMatrix(gm, 0, k);
}

/*
# `lastColumns`

Return the last `k` columns of `gm` with no copying.

**Warning:** See the warning on `middleColumns`.
*/

GretlMatrix lastColumns(GretlMatrix gm, int k) {
	return GretlMatrix(gm, gm.cols-k, gm.cols);
}

/*
# sort

Sort the rows of matrix `m` according to column `col`. By default, the
sort is done in ascending order, but you can change `direction` to
`"dec"` to sort in decreasing order.

*Example:*

```
m.sort(2); // Sort in ascending order using the third column
m.sort!"dec"(4) // Sort in decreasing order using the fifth column
```

*/

private struct Observation {
  double value;
  int obs;
}

DoubleMatrix sort(string direction="inc")(GretlMatrix m, int col) {
  // Put all elements and row numbers in an array
  Observation[] temp;
  foreach(ii; 0..m.rows) { 
    temp ~= Observation(m[ii,col], ii); 
  }

  // Sort the array according to data
  static if (direction == "inc") {
    bool com(Observation x, Observation y) { return x.value < y.value; }
  } else static if (direction == "dec") {
    bool com(Observation x, Observation y) { return x.value > y.value; }
  } else {
		static assert(false, `The only options for sort direction are "inc" and "dec".`);
	}
  auto output = std.algorithm.sort!(com)(temp);

  // Now put the results into a new matrix
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows) { Row(result, ii) = Row(m, output[ii].obs); }
  return result;
}

int[] order(bool inc=true)(GretlMatrix m, int col) {
  Observation[] temp;
  foreach(ii; 0..m.rows) { 
    temp ~= Observation(m[ii,col], ii);
  }

  // Sort the array according to data
  static if (inc) {
    bool com(Observation x, Observation y) { return x.value < y.value; }
  } else {
    bool com(Observation x, Observation y) { return x.value > y.value; }
  }
  auto output = std.algorithm.sort!(com)(temp);
  int[] result;
  foreach(val; output) {
                result ~= val.obs;
        }
  return result;
}

DoubleMatrix stackColumns(GretlMatrix m) {
  auto result = DoubleMatrix(m.rows*m.cols, 1);
  foreach(ii; 0..result.rows) {
    result[ii,0] = m.ptr[ii];
  }
  return result;
}

DoubleMatrix stackRows(GretlMatrix m) {
  return stackColumns(m.transpose);
}

T cbind(T)(T[] mats) {
  int rows = mats[0].rows;
  int cols = 0;
  foreach(mat; mats) {
    enforce(mat.rows == rows, "Number of rows is not the same for all arguments");
    cols += mat.cols;
  }
  auto result = T(rows, cols);
  int startColumn = 0;
  foreach(mat; mats) {
    auto temp = result[0..rows, startColumn..startColumn+mat.cols];
    temp = Submatrix(mat);
    startColumn += mat.cols;
  }
  return result;
}

DoubleMatrix cbind(double[][] arrays) {
	DoubleMatrix[] tmp;
	foreach(arr; arrays) {
		tmp ~= DoubleMatrix(arr);
	}
	return cbind!DoubleMatrix(tmp);
}

T rbind(T)(T[] mats) {
  int cols = mats[0].cols;
  int rows = 0;
  foreach(mat; mats) {
    enforce(mat.cols == cols, "Number of cols is not the same for all arguments");
    rows += mat.rows;
  }
  auto result = T(rows, cols);
  int startRow = 0;
  foreach(mat; mats) {
    result[startRow..startRow+mat.rows, 0..mat.cols] = mat;
    startRow += mat.rows;
  }
  return result;
}

DoubleMatrix rbind(double[][] arrays) {
	DoubleMatrix[] tmp;
	foreach(arr; arrays) {
		tmp ~= DoubleMatrix(arr);
	}
	return rbind!DoubleMatrix(tmp);
}

DoubleMatrix each(alias f)(GretlMatrix y, GretlMatrix x) {
  auto result = DoubleMatrix(y.rows, y.cols);
  foreach(ii; 0..y.rows*y.cols) {
    result.ptr[ii] = f(y.ptr[ii], x.ptr[ii]);
  }
  return result;
}

DoubleMatrix each(alias f)(GretlMatrix m) {
  auto result = DoubleMatrix(m.rows, m.cols);
  foreach(ii; 0..m.rows*m.cols) {
    result.ptr[ii] = f(m.ptr[ii]);
  }
  return result;
}

// Fill the given observations with val
void fillObs(GretlMatrix m, int[] rows, int[] cols, double val) {
  enforce(rows.length == cols.length, "Need the same number of row and column observations");
  foreach(row, col; lockstep(rows, cols)) {
    m[row, col] = val;
  }
}

DoubleMatrix selectColumns(GretlMatrix x, int[] columns) {
	auto result = DoubleMatrix(x.rows, columns.sum);
	int fillPointer = 0;
	foreach(col, int flag; columns) {
		if (flag == 1) {
			Col(result, fillPointer) = Col(x, col.to!int);
			fillPointer += 1;
		}
	}
	return result;
}

// For use inside structs, for which opAssign doesn't work
void replace(ref DoubleMatrix m, DoubleMatrix newmat) {
	m.data = newmat.data;
	m.rows = newmat.rows;
	m.cols = newmat.cols;
	//m.data = new double[newmat.rows*newmat.cols];
	//m.data[] = newmat.data[];
}
