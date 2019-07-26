/* Basically the same example as logit.d, but I want to show how to do it. */
import gretl.base, gretl.random;

void main() {
	randInit();
	setSeed(200); // To compare with logit
	
	// Dataset with 4 variables, column 0 is the intercept
	// Columns 1 and 2 are
	// Initialize the dataset by calling create_new_dataset
	// Column 0 is initialized to 1 for the intercept
	// Use a scope statement to free the memory allocated by Gretl
	// Will do that when done, but can put it here so we don't forget
	Dataset * ds = create_new_dataset(4, 100, 0);
	scope(exit) { destroy_dataset(ds); }
	
	// Dataset has opIndex defined, so you can treat it like a matrix
	// to avoid bugs, admittedly with suboptimal notation
	//~ foreach(ii; 0..100) {
		//~ ds.Z[1][ii] = rnorm();
		//~ ds.Z[2][ii] = rnorm();
	//~ }
	foreach(ii; 0..100) {
		(*ds)[ii,1] = rnorm();
		(*ds)[ii,2] = rnorm();
	}
	// Haven't yet implemented multidimensional slicing
	// Can still take advantage of D's built-in slicing
	ds.Z[3][0..50] = 1.0;
	ds.Z[3][50..100] = 0.0;

	// Specify the model using an array of ints
	// 4 elements, third variable is the lhs variable, 0 (intercept),
	// 1, and 2 are the regressors
	int[] spec = [4, 3, 0, 1, 2];
	
	// Initialize prn by calling gretl_print_new
	GretlPrn * prn = gretl_print_new(PrintType.stdout, null);
	scope(exit) { gretl_print_destroy(prn); }
	
	// It's easy to mess this up
	// Initialize logitModel by calling gretl_model_new
	Model * probitModel = gretl_model_new();
	scope(exit) { gretl_model_free(probitModel); }
	*probitModel = binary_probit(spec.ptr, ds, gretlopt.none, prn);
	
	// Handy print functions for all of the models
	printmodel(probitModel, ds, gretlopt.none, prn);
}
