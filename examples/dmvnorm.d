/* Test the dmvnorm results in here against the mvtnorm
 * package for R. Exact same results in each case. It's
 * possible that edge cases will be handled differently,
 * but this is good enough for me.*/
module mvnormtest;

import gretl.matrix, gretl.vector, gretl.random;
import std.stdio;

void main() {
	auto x = DoubleVector([0.1, 0.2, 0.3]);
	auto mu = DoubleVector([0.1, 0.2, 0.3]);
	auto V = DoubleMatrix(3,3);
	Row(V,0) = [1.0, 0.0, 0.0];
	Row(V,1) = [0.0, 1.0, 0.0];
	Row(V,2) = [0.0, 0.0, 1.0];
	inv(V).print("Inverse:");
	writeln(dmvnorm(x, mu, V));
	writeln(dmvnorm(x, mu, V, true));

	auto mu2 = DoubleVector([1.0, 4.2, 0.3]);
	writeln(dmvnorm(x, mu2, V));
	writeln(dmvnorm(x, mu2, V, true));
}

/* R test code
 * mvtnorm version 1.1-1
library(mvtnorm)
v <- matrix(c(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0), nrow=3)
print(dmvnorm(c(0.1, 0.2, 0.3), c(0.1, 0.2, 0.3), v, log=FALSE))
print(dmvnorm(c(0.1, 0.2, 0.3), c(1.0, 4.2, 0.3), v, log=FALSE))
print(dmvnorm(c(0.1, 0.2, 0.3), c(0.1, 0.2, 0.3), v, log=TRUE))
print(dmvnorm(c(0.1, 0.2, 0.3), c(1.0, 4.2, 0.3), v, log=TRUE))
*/
