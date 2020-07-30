import std.stdio;
import gretl.matrix;
import std.stdio;

void main() {
	auto m = DoubleMatrix(3,3);

	/* Fill using Row structs */
	Row(m, 0) = [1.1, 2.2, 3.3];
	Row(m, 1) = [4.4, 5.5, 6.6];
	Row(m, 2) = [7.7, 8.8, 9.9];
	
	/* Print */
	gretl.matfunctions.print(m, "(3x3) matrix m");

	/* Get a submatrix and print it out */
	m[0..2, 0..2].print("Upper left corner of m");
	m[1..3, 1..2].print("Some elements");
	
	/* Do matrix multiplication with Submatrix structs */
	DoubleMatrix m2 = m[0..2, 1..3] * m[1..3, 0..2];
	m2.print("Submatrix product");
	
	/* Convert Submatrix to double[] */
	writeln(m[_all,0].array);
	writeln(m[1..3, 2].array);
	writeln(m[1, 0..3].array);
	
	/* Slicing of a Submatrix */
	Submatrix m3 = m[_all, 0];
	writeln(m3[0..2]);
	Submatrix m4 = m[1, _all];
	writeln(m4[1..3]);
	
	/* Slicing of Row and Col */
	writeln(Row(m,1)[1..3]);
	writeln(Col(m,0)[0..2]);
}

