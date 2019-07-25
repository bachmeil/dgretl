import gretl.matrix;
import std.stdio;

void main() {
    auto m = DoubleMatrix(3,3);
    m.print("Introducing our new matrix m");
    
    /* Fill using Row structs */
    Row(m, 0) = [1.1, 2.2, 3.3];
    Row(m, 1) = [4.4, 5.5, 6.6];
    Row(m, 2) = [7.7, 8.8, 9.9];
    
    /* Print */
    m.print("(3x3) matrix m");
}
