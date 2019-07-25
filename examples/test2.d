import std.stdio;
import gretl.random;

void main() {
	// Have to initialize the random number generator
	randInit();
	
	// Can work with the seed
	writeln("The seed is ", getSeed());
	setSeed(250);
	writeln("The seed is ", getSeed());
	
	// Generate integers between 0 and 2
	foreach(_; 0..6) {
		writeln(genInteger(3));
	}	
	
	// Generate uniform random numbers between 0 and 1
	foreach(_; 0..6) {
		writeln(runif());
	}
	
	// Alternatively, generate U[0,1] vector
	runif(10).print("U[0,1] vector");
	
	// Generate U[10, 30] vector
	runif(10, 10.0, 30.0).print("U[10, 30] vector");
	
	// Generate N(0,1) variables
	foreach(_; 0..6) {
		writeln(rnorm());
	}
	
	// Vector of N(0,1)
	rnorm(10).print("N(0,1) vector");
	
	// Vector of N(1.0, 5.0)
	rnorm(10, 1.0, 5.0).print("N(1.0, 5.0) vector");
}
