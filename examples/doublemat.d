import gretl.base, gretl.matrix;

void main() {
  /* Matrix construction */
  auto dm1 = DoubleMatrix(5); // Default is one column
  auto dm2 = DoubleMatrix(5,2); // (5x2)
