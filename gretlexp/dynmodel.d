/* Experimental support for time series data
 * This was in a previous version, but it needs to be totally reworked.
 * 
 * November 30, 2017
 * 
 * For now, TS objects will store data in a DoubleMatrix. That means it's
 * managed by the D garbage collector. The cost is that any data you pass
 * from R to a TS object requires copying.
 * 
 * Initially, I considered having this be a reference to the data, or
 * at least having the option to manually manage the memory. That doesn't
 * work, because every time you do an arithmetic operation, take the difference
 * or log of a variable, and so on, you create a new TS struct. That's
 * beyond the target audience of this library, but it wouldn't take much
 * to do that if it's something you want for yourself.
 * 
 * This library produces fast code relative to R or Matlab, but it is
 * going to be slower than an optimized C library, because the priority
 * is ease of use.
 */
module gretlexp.dynmodel;
 
import std.conv, std.exception, std.math, std.stdio, std.utf;
 
import gretl.base, gretl.matrix, gretl.reg;
version(r) {
	import embedr.r;
}

enum structure {
	cross_section,
	time_series,
	stacked_time_series,
	stacked_cross_section,
	panel_unknown,
	special_time_series,
	unknown
};

struct Index {
	int ii;
	
	this(int val) {
		ii = val;
	}
	
	this(ulong val) {
		ii = val.to!int;
	}
}

private ulong[2] toDate(ulong x, ulong f) {
	ulong y = (x-1)/f;
	ulong m = x-f*y;
	return [y, m];
}

private ulong toInt(ulong[2] d, ulong f) {
	return f*d[0] + d[1];
}

// Create a Dataset for use in Gretl
// Contains a reference to the data
// Intended to be immediately destroyed after use
Dataset asDataset(TS x, double** z) {
	Dataset result;
	result.v = 1;
	result.n = x.length.to!int;
	result.Z = z;
	dataset_set_time_series(&result, x.frequency.to!int, x.start[0].to!int, x.start[1].to!int);
	return result;
}

// Needs manual freeing and also copies
Dataset * mallocDataset(TS x) {
	Dataset * result = create_new_dataset(1, x.length.to!int, 0);
	dataset_set_time_series(result, x.frequency.to!int, x.start[0].to!int, x.start[1].to!int);
	foreach(ii; 0..x.length) {
		(*result)[ii.to!int,0] = x.mat.data.ptr[ii];
	}
	return result;
}

/**
 * Holds the output of a call into Gretl to estimate
 * an ARMA model.
 */
struct ArmaFit {
	int p;
	int q;
	double[] coef;
	double[] sderr;
	DoubleMatrix v;
	int nobs;
	double ess;
	double sigma;
	double lnL;
	double ybar;
	double sdy;
	double[string] criteria;
	
	this(Model * m, int _p, int _q) {
		p = _p;
		q = _q;
		foreach(ii; 0..m.ncoeff) {
			coef ~= m.coeff[ii];
			sderr ~= m.sderr[ii];
		}
		nobs = m.nobs;
		ess = m.ess;
		sigma = m.sigma;
		lnL = m.lnL;
		ybar = m.ybar;
		sdy = m.sdy;
		addCriteria();
		addV(m);
	}		
	
	void addCriteria() {
		criteria["aic"] = -2.0 * lnL + 2.0 * coef.length;
		criteria["sic"] = -2.0 * lnL + std.math.log(nobs) * coef.length;
		criteria["hqc"] = -2.0 * lnL + 2  * std.math.log(std.math.log(nobs)) * coef.length;
	}
	
	void addV(Model * m) {
		v = DoubleMatrix(coef.length.to!int, coef.length.to!int);
		foreach(col; 0..v.cols) {
			foreach(row; 0..v.rows) {
				v[row, col] = gretl_model_get_vcv_element(m, row, col, coef.length.to!int);
			}
		}
	}
	
	void print() {
		writeln("Estimated ARMA(", p, ",", q, ") Model");
		writeln("Intercept:     ", coef[0]);
		writeln("sd(Intercept): ", sderr[0]);
		writeln("");
		if (p > 0) {
			writeln("AR coefficients:     ", coef[1..p+1]);
			writeln("sd(AR coefficients): ", sderr[1..p+1]);
		}
		if (q > 0) {
			writeln("MA coefficients:     ", coef[p+1..$]);
			writeln("sd(MA coefficients): ", sderr[p+1..$]);
		}
		gretl.matrix.print(v, "Residual covariance matrix:");
		writeln("ESS: ", ess);
		writeln("sigma: ", sigma);
		writeln("log L: ", lnL);
		writeln("Information criteria for the estimated model:");
		writeln(criteria);
	}
}

// Copy the data
//~ ArmaFit armaGretl(TS x, ulong p, ulong q) {
	//~ int[] list = [5, p.to!int, 0, q.to!int, gretlListsep, 0];
	//~ Dataset * d = create_new_dataset(1, x.length.to!int, 0);
	//~ foreach(ii; 0..x.length.to!int) {
		//~ d.Z[0][ii] = x.mat[ii,0];
	//~ }
	//~ dataset_set_time_series(d, x.frequency.to!int, x.start[0].to!int, x.start[1].to!int);
	//~ GretlPrn * prn = gretl_print_new(PrintType.stdout, null);
	//~ Model *m;
	//~ m = gretl_model_new();
	//~ *m = arma(list.ptr, null, d, gretlopt.none.to!int, prn);
	//~ auto result = ArmaFit(m, p.to!int, q.to!int);
	//~ gretl_model_free(m);
	//~ destroy_dataset(d);
	//~ return result;
//~ }

// This works for some reason without copying
ArmaFit armaGretl(TS x, ulong p, ulong q) {
	int[] list = [5, p.to!int, 0, q.to!int, gretlListsep, 0];
	Dataset * d = create_new_dataset(1, 1, 0);
	double** tmp = d.Z;
	double*[1] z;
	z[0] = x.mat.data.ptr;
	d.Z = z.ptr;
	d.n = x.length.to!int;
	d.t2 = d.n-1;
	dataset_set_time_series(d, x.frequency.to!int, x.start[0].to!int, x.start[1].to!int);
	GretlPrn * prn = gretl_print_new(PrintType.stdout, null);
	Model *m;
	m = gretl_model_new();
	*m = arma(list.ptr, null, d, gretlopt.none.to!int, prn);
	auto result = ArmaFit(m, p.to!int, q.to!int);
	d.Z = tmp;
	gretl_model_free(m);
	destroy_dataset(d);
	return result;
}

// This always fails for some reason
//~ ArmaFit armaGretl(TS x, ulong p, ulong q) {
	//~ int[] list = [5, p.to!int, 0, q.to!int, gretlListsep, 0];
	//~ Dataset * d = new Dataset;
	//~ writeln("d.Z");
	//~ writeln(d.Z);
	//~ double*[1] z;
	//~ z[0] = x.mat.data.ptr;
	//~ d.Z = z.ptr;
	//~ writeln(d.Z);
	//~ writeln(d.Z[0][5]);
	//~ d.v = 1;
	//~ d.n = x.length.to!int;
	//~ d.t2 = x.length.to!int-1;
	//~ d.structure = 0;
	//~ d.pd = 1;
	//~ d.sd0 = 1.0;
	//~ d.stobs[] = '\0';
	//~ d.stobs[0] = '1';
	//~ d.endobs[] = '\0';
	//~ d.endobs[0] = '1';
	//~ d.markers ='\0';
	//~ d.modflag = '\0';
	//~ d.panel_sd0 = 0.0;
	//~ dataset_set_time_series(d, x.frequency.to!int, x.start[0].to!int, x.start[1].to!int);

	//~ GretlPrn * prn = gretl_print_new(PrintType.stdout, null);
	//~ writeln("before");
	//~ // This is where it fails
	//~ // d is exactly the same as in the other function
	//~ Model * m;
	//~ m = gretl_model_new();
	//~ *m = arma(list.ptr, null, d, gretlopt.none.to!int, prn);
	//~ writeln("after");
	//~ auto result = ArmaFit(m, p.to!int, q.to!int);
	//~ gretl_model_free(m);
	//~ return result;
//~ }

/* For now, we can use the second time element as zero for annual data */
struct TS {
	DoubleMatrix mat;
	ulong offset;
	ulong frequency;

	// Copies
	this(GretlMatrix m, ulong[2] s, ulong f) {
		enforce(m.cols == 1, "Cannot convert to a TS object unless there is only one column. You may want to use a TS[] array instead.");
		mat = DoubleMatrix(m.rows, 1);
		foreach(ii; 0..m.rows) {
			mat.data[ii] = m.ptr[ii];
		}
		offset = s.toInt(f);
		frequency = f;
	}
	
	// Copies
	this(DoubleVector v, ulong[2] s, ulong f) {
		mat = DoubleMatrix(v.length, 1);
		foreach(ii; 0..v.length) {
			mat.data[ii] = v[ii];
		}
		offset = s.toInt(f);
		frequency = f;
	}

	this(ulong[2] s, ulong[2] e, ulong f) {
		ulong obs = e.toInt(f) - s.toInt(f) + 1;
		mat = DoubleMatrix(obs.to!int, 1);
		offset = s.toInt(f);
		frequency = f;
	}
	
	this(TS x, ulong[2] s, ulong[2] e) {
		enforce(s.toInt(x.frequency) >= x.offset, "Start date is before the first observation");
		enforce(e.toInt(x.frequency) <= x.end.toInt(x.frequency), "End data is after the last observation");
		mat.data = x.mat.data[s.toInt(x.frequency) - x.offset..$];
		offset = s.toInt(x.frequency);
		frequency = x.frequency;
	}
	
	this(double[] v, ulong[2] s, ulong f) {
		mat = DoubleMatrix(v.length.to!int, 1);
		foreach(ii; 0..v.length.to!int) {
			mat.data[ii] = v[ii];
		}
		offset = toInt(s, f);
		frequency = f;
	}

	this(double x, ulong[2] s, ulong f) {
		this([x], s, f);
	}
	
	double * ptr() {
		return mat.data.ptr;
	}
	
	/* Want to allow access by index
	 * But also want to be explicit about it to prevent bugs */
	double opIndex(Index ind) {
		enforce(ind.ii >= 0, "Cannot have a negative index on a TS struct");
		enforce(ind.ii < mat.rows, "Index on TS struct too large");
		return mat.data[ind.ii];
	}
	
	/* This function does bounds checking and returns the index number
	 * associated with the given date */
	ulong getIndex(ulong y, ulong m) {
		ulong v = y*frequency + m;
		ulong ind = v - offset;
		enforce(ind >= 0, "Trying to access elements prior to the start of the dataset");
		enforce(ind < mat.rows, "Trying to access elements after the end of the dataset");
		return ind;
	}
	
	double opIndex(ulong y, ulong m) {
		return mat.data[getIndex(y, m)];
	}
	
	// Pass date as an int rather than ulong[2]
	// This is NOT the same as passing an index number!
	double intIndex(ulong d) {
		return mat.data[d-offset];
	}
	
	void opIndexAssign(double val, ulong y, ulong m) {
		mat.data[getIndex(y, m)] = val;
	}
	
	// Pass date as an int rather than ulong[2]
	// This is NOT the same as passing an index number!
	void intIndexAssign(double val, ulong d) {
		mat.data[d-offset] = val;
	}

	void opIndexAssign(double val, Index ind) {
		enforce(ind.ii >= 0, "Cannot have a negative index on a TS struct");
		enforce(ind.ii < mat.rows, "Index on TS struct too large");
		mat.data[ind.ii] = val;
	}
	
	TS opBinary(string op)(TS y) {
		static if(op == "+") {
			return addition(this, y);
		}
		static if(op == "-") {
			return subtraction(this, y);
		}
		static if(op == "*") {
			return multiplication(this, y);
		}
		static if(op == "/") {
			return division(this, y);
		}
	}
	
	TS opBinary(string op)(double a) {
		static if(op == "+") {
			return addition(this, a);
		}
		static if(op == "-") {
			return subtraction(this, a);
		}
		static if(op == "*") {
			return multiplication(this, a);
		}
		static if(op == "/") {
			return division(this, a);
		}
	}
	
	TS opBinaryRight(string op)(double a) {
		static if(op == "+") {
			return addition(this, a);
		}
		static if(op == "-") {
			return subtraction(a, this);
		}
		static if(op == "*") {
			return multiplication(this, a);
		}
		static if(op == "/") {
			return division(a, this);
		}
	}

	ulong[2] start() {
		return toDate(offset, frequency);
	}
	
	ulong intStart() {
		return offset;
	}
	
	ulong[2] end() {
		return toDate(offset+mat.rows-1, frequency);
	}
	
	ulong intEnd() {
		return offset+mat.rows-1;
	}
	
	ulong length() {
		return mat.rows.to!ulong;
	}
	
	TS opSlice(ulong[2] s, ulong[2] e) {
		enforce(s.toInt(frequency) >= offset, "first index is before the start of the time series");
		enforce(e.toInt(frequency) <= intEnd(), "second index is after the end of the time series");
		TS result;
		result.mat = mat;
		result.mat.data = mat.data[s.toInt(frequency)-offset..e.toInt(frequency)-offset+1];
		result.mat.rows = to!int(e.toInt(frequency) - s.toInt(frequency) + 1);
		result.mat.cols = 1;
		result.frequency = frequency;
		result.offset = s.toInt(frequency);
		return result;
	}
	
	DoubleVector vec() {
		return DoubleVector(mat);
	}
	
	//~ void print() {
		//~ writeln("Start date: ", start());
		//~ writeln("End date: ", end());
		//~ writeln("Frequency: ", frequency);
	//~ }
}

void print(TS x) {
	foreach(obs; 0..x.length) {
		ulong[2] date = toDate(obs+x.offset, x.frequency);
		string space = " 0";
		if (date[1] > 9) { space = " "; }
		writeln(date[0], space, date[1], " ", x[Index(obs.to!int)]);
	}
}

ulong[2] maxDate(ulong[2] date1, ulong[2] date2) {
	if (date1[0] > date2[0]) {
		return date1;
	} else if (date2[0] > date1[0]) {
		return date2;
	} else {
		if (date1[1] > date2[1]) {
			return date1;
		} else {
			return date2;
		}
	}
}

ulong[2] lastStart(MTS vars) {
	enforce(vars.length > 0, "Can't get the latest start date of an MTS struct with no elements");
	ulong[2] result = vars[0].start;
	if (vars.length > 1) {
		foreach(var; vars[1..$]) {
			result = maxDate(result, var.start);
		}
	}
	return result;
}

ulong[2] lastStart(TS y, MTS x) {
	x ~= y;
	return lastStart(x);
}

ulong[2] minDate(ulong[2] date1, ulong[2] date2) {
	if (date1[0] < date2[0]) {
		return date1;
	} else if (date2[0] < date1[0]) {
		return date2;
	} else {
		if (date1[1] < date2[1]) {
			return date1;
		} else {
			return date2;
		}
	}
}

ulong[2] firstEnd(MTS vars) {
	enforce(vars.length > 0, "Can't get the earliest end date of an MTS struct with no elements");
	ulong[2] result = vars[0].end;
	if (vars.length > 1) {
		foreach(var; vars[1..$]) {
			result = minDate(result, var.end);
		}
	}
	return result;
}

ulong[2] firstEnd(TS y, MTS x) {
	x ~= y;
	return firstEnd(x);
}

TS addition(TS y, TS x) {
	enforce(y.frequency == x.frequency, "Cannot add TS structs of different frequency");
	auto result = TS(maxDate(y.start, x.start), minDate(y.end, x.end), y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(y.intIndex(ii) + x.intIndex(ii), ii);
	}
	return result;
}

TS addition(TS y, double a) {
	auto result = TS(y.start, y.end, y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(y.intIndex(ii) + a, ii);
	}
	return result;
}

TS subtraction(TS y, TS x) {
	enforce(y.frequency == x.frequency, "Cannot add TS structs of different frequency");
	auto result = TS(maxDate(y.start, x.start), minDate(y.end, x.end), y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(y.intIndex(ii) - x.intIndex(ii), ii);
	}
	return result;
}

TS subtraction(TS y, double a) {
	auto result = TS(y.start, y.end, y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(y.intIndex(ii) - a, ii);
	}
	return result;
}

TS subtraction(double a, TS y) {
	auto result = TS(y.start, y.end, y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(a - y.intIndex(ii), ii);
	}
	return result;
}

TS multiplication(TS y, TS x) {
	enforce(y.frequency == x.frequency, "Cannot add TS structs of different frequency");
	auto result = TS(maxDate(y.start, x.start), minDate(y.end, x.end), y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(y.intIndex(ii) * x.intIndex(ii), ii);
	}
	return result;
}

TS multiplication(TS y, double a) {
	auto result = TS(y.start, y.end, y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(y.intIndex(ii) * a, ii);
	}
	return result;
}

TS division(TS y, TS x) {
	enforce(y.frequency == x.frequency, "Cannot add TS structs of different frequency");
	auto result = TS(maxDate(y.start, x.start), minDate(y.end, x.end), y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(y.intIndex(ii) / x.intIndex(ii), ii);
	}
	return result;
}

TS division(TS y, double a) {
	auto result = TS(y.start, y.end, y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(y.intIndex(ii) / a, ii);
	}
	return result;
}

TS division(double a, TS y) {
	auto result = TS(y.start, y.end, y.frequency);
	foreach(ii; result.intStart..result.intEnd+1) {
		result.intIndexAssign(a / y.intIndex(ii), ii);
	}
	return result;
}

TS lag(TS x, ulong k) {
	TS result;
	result.mat = x.mat;
	result.offset = x.offset+k;
	result.frequency = x.frequency;
	return result;
}

TS diff(TS x, ulong k=1) {
	return x - lag(x, k);
}

TS trend(TS x, ulong k=1) {
	TS result = TS(x.start, x.end, x.frequency);
	foreach(ii; 0..result.length) {
		result[Index(ii)] = (ii+1)^^k;
	}
	return result;
}

TS log(TS x) {
	TS result = TS(x.start, x.end, x.frequency);
	foreach(ii; 0..result.length) {
		result[Index(ii)] = std.math.log(x[Index(ii)]);
	}
	return result;
}

long nobs(ulong[2] s, ulong[2] e, ulong f) {
	return e.toInt(f) - s.toInt(f);
}

struct ModelVar {
	DoubleVector y;
	DoubleMatrix x;
}

struct TSModel {
	TS lhs;
	MTS rhs;
	
	this(TS y, TS[] x) {
		lhs = y;
		rhs = MTS(x);
	}
	
	this(TS y, MTS x) {
		lhs = y;
		rhs = x;
	}
	
	ModelVar regressors() {
		ModelVar result;
		ulong[2] start = lastStart(lhs, rhs);
		ulong[2] end = firstEnd(lhs, rhs);
		TS depvar = lhs[start..end];
		result.y = depvar.vec;
		result.x = DoubleMatrix(depvar.length.to!int, rhs.length.to!int);
		foreach(col; 0..result.x.cols) {
			TS x = rhs[col];
			foreach(row; 0..result.x.rows) {
				result.x[row, col] = x.mat[row,0];
			}
		}
		return result;
	}
	
	//void addIntercept() {}
	//Find the minimum start date and maximum end date
	//Create a ts of ones for those dates
}

OlsFit tsreg(TSModel mod) {
	ModelVar yx = mod.regressors;
	return lm(yx.y, yx.x);
}

struct MTS {
	TS[] data;
	
	alias data this;
	
	this(TS x) {
		data = [x];
	}
	
	this(TS[] x) {
		data = x;
	}
	
	void put(TS x) {
		data ~= x;
	}
	
	void put(TS[] xs) {
		foreach(x; xs) {
			data ~= x;
		}
	}
	
	void opOpAssign(string op)(MTS x) {
		static if(op == "~") {
			this.data ~= x.data;
		}
	}
	
	void opOpAssign(string op)(TS x) {
		static if(op == "~") {
			this.data ~= x;
		}
	}
	
	//DoubleMatrix array() {}
	
	/* When you know with certainty the start and end dates */
	//~ DoubleMatrix array(ulong[2] s, ulong[2] e) {
		//~ auto result = DoubleMatrix(nobs(s, e, data[0].frequency).to!int, data.length.to!int);
		//~ foreach(ii, tmp; data) {
			//~ Col(result, ii.to!int) = tmp;
		//~ }
		//~ return result;
	//~ }
}

MTS lags(TS x, long[] ks) {
	MTS result;
	foreach(k; ks) {
		result.put(lag(x,k));
	}
	return result;
}

MTS lags(TS x, long k0, long k1) {
	MTS result;
	foreach(k; k0..k1+1) {
		result.put(lag(x,k));
	}
	return result;
}

MTS lags(TS x, long k) {
	MTS result;
	result.put(lag(x,k));
	return result;
}

/* Annual time series should be treated differently because the
 * date is identified by a single number, while at higher frequencies
 * you need two numbers. Alternatively, I could convert everything to
 * floating point, but that's a mess.
 */
//~ struct AnnualTS {
	//~ TS x;
	//~ alias x this;

	//~ this(GretlMatrix m, ulong s) {
		//~ enforce(m.cols == 1, "Cannot create AnnualTS struct if there is more than one column. You may want to use an ATS[] array instead.");
		//~ x.data.rows = m.rows;
		//~ x.data.cols = 1;
		//~ //x.data.ptr = m.ptr;
		//~ x.offset = s;
		//~ x.frequency = 1;
	//~ }
	
	//~ ulong getIndex(ulong y) {
		//~ ulong ind = y - x.offset;
		//~ enforce(ind >= 0, "Trying to access elements prior to the start of the dataset");
		//~ enforce(ind < x.data.rows, "Trying to access elements after the end of the dataset");
		//~ return ind;
	//~ }

	//~ double opIndex(Index ind) {
		//~ return x.opIndex(ind); 
	//~ }
	
	//~ double opIndex(ulong y) {
		//~ return x.data.ptr[getIndex(y)];
	//~ }
	
	//~ void opIndexAssign(double val, ulong y) {
		//~ x.data.ptr[getIndex(y)] = val;
	//~ }
	
	//~ ulong start() {
		//~ return x.offset;
	//~ }

	//~ ulong end() {
		//~ return x.offset+x.data.rows-1;
	//~ }
	
	//~ // length is provided by TS
//~ }
