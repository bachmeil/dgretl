module gretl.base;

import gretl.matrix;
import std.utf;
version(r) {
	import embedr.r;
}

immutable int gretlListsep = -100;
private static immutable int maxline = 16384;
private static immutable int maxlabel = 128;
private static immutable int maxlen = 512;
private static immutable int maxdisp = 32;
private static immutable int vnamelen = 32;
private static immutable int obslen = 16;

private struct ARInfo {
	int * arlist;
	double * rho;
	double * sderr;
}

// GretlMatrixMod
enum matmod {
    none = 0,
    transpose,
    square,
    cumulate,
    decrement
}

// GretlMatrixStructure
enum matstructure {
    square = 1,
    lowertri,
    uppertri,
    symmetric,
    diagonal,
    identity,
    scalar
}

enum gretlopt {
	none = 0,
	a = 1 << 0,
	b = 1 << 1,
	c = 1 << 2,
	d = 1 << 3,
	e = 1 << 4,
	f = 1 << 5,
	g = 1 << 6,
	h = 1 << 7,
	i = 1 << 8,
	j = 1 << 9,
	k = 1 << 10,
	l = 1 << 11,
	m = 1 << 12,
	n = 1 << 13,
	o = 1 << 14,
	p = 1 << 15,
	q = 1 << 16,
	r = 1 << 17,
	s = 1 << 18,
	t = 1 << 19,
	u = 1 << 20,
	v = 1 << 21,
	w = 1 << 22,
	x = 1 << 23,
	z = 1 << 24,
	y = 1 << 25,
	unset = 1 << 30
}

enum GretlFileType {
    xml,       /* gretl XML data file (.gdt) */
    binary,    /* zip file with binary component (.gdtb) */
    csv,            /* comma-separated or other plain text data */
    ocatve,         /* GNU octave ascii data file */
    gnumeric,       /* gnumeric workbook data */
    xls,            /* MS Excel spreadsheet data */
    xlsx,           /* MS Office Open XML spreadsheet data */
    ods,            /* Open Document Spreadsheet data */
    wf1,            /* Eviews workfile data */
    dta,            /* Stata .dta data */
    sav,            /* SPSS .sav data */
    sas,            /* SAS xport data file */
    jmulti,         /* JMulTi data file */
    max,       /* -- place marker -- */
    script,         /* file containing gretl commands */
    session,        /* zipped session file */
    db,      /* gretl database */
    dbwww,  /* gretl database, accessed via internet */
    rats,        /* RATS 4.0 database */
    pcgive,      /* PcGive bn7/in7 pair */
    odbc,           /* Open DataBase Connectivity */
    other    /* none of the above */
}

enum ModelSelCriteria {
	aic,
	bic,
	hqc,
	max,
}

enum PrintType {
    stdout,
    stderr,
    file,
    buffer,
    tempfile,
    stream
}

private struct Sample {
	int t1;
	int t2;
	uint rseed;
}

private struct ModelTest{}
private struct ModelDataItem{}
private struct VarInfo{}
struct GretlPrn{}

struct Dataset {
  int v;              /* number of variables */
  int n;              /* number of observations */
  int pd;             /* periodicity or frequency of data */
  int structure;      /* time series, cross section or whatever */
  double sd0;         /* float representation of stobs */
  int t1, t2;         /* start and end of current sample */
  char[obslen] stobs;  /* string representation of starting obs (date) */
  char[obslen] endobs; /* string representation of ending obs */
  double **Z;         /* data array */
  char **varname;     /* array of names of variables */
  VarInfo **varinfo;  /* array of specific info on vars */
  char markers;       /* NO_MARKERS (0), REGULAR MARKERS or DAILY_DATE_STRINGS */
  version(old) {}
  else {
		char modflag;       /* binary flag for dataset modified or not */
	}
  char **S;           /* to hold observation markers */
  char *descrip;      /* to hold info on data sources etc. */
  char *submask;      /* subsampling mask */
  char *restriction;  /* record of sub-sampling restriction */
  char *padmask;      /* record of padding to re-balance panel data */
  version(old) {}
  else {
		uint rseed; /* resampling seed */
  }
  int auxiliary;      /* = 0 for regular dataset, 1 for aux dataset */
  char *pangrps;      /* panel-only: name of series holding group names */
  int panel_pd;       /* panel-only: panel time-series frequency */
  double panel_sd0;   /* panel-only: time-series start */
  
  double opIndex(int row, int col) {
		return Z[col][row];
	}
	
	void opIndexAssign(double val, int row, int col) {
		Z[col][row] = val;
	}
}

struct Model {
  int ID;                      /* ID number for model */
  int refcount;                /* for saving/deleting */
  int ci;                      /* "command index" -- estimation method */
  gretlopt opt;                /* record of options */
  int t1, t2, nobs;            /* starting observation, ending
                                  observation, and number of obs */
  char *submask;               /* keep track of sub-sample in force
                                  when model was estimated */
  char *missmask;              /* missing observations mask */
  Sample smpl;                 /* numeric start and end of current sample
                                    when model was estimated */
  version(old) {}
  else {
    int full_n;                  /* full length of dataset on estimation */
	}
  int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                  freedom in numerator and denominator */
  int *list;                   /* list of variables by ID number */
  int ifc;                     /* = 1 if the equation includes a constant,
	    else = 0 */
  int nwt;                     /* ID number of the weight variable (WLS) */
  int aux;                     /* code representing the sort of
				    auxiliary regression this is (or not) */
  double *coeff;               /* array of coefficient estimates */
  double *sderr;               /* array of estimated std. errors */
  double *uhat;                /* regression residuals */
  double *yhat;                /* fitted values from regression */
  double *xpx;                 /* X'X matrix, in packed form */
  double *vcv;                 /* VCV matrix for coefficient estimates */
  double ess, tss;             /* Error and Total Sums of Squares */
  double sigma;                /* Standard error of regression */
  double rsq, adjrsq;          /* Unadjusted and adjusted R^2 */     
  double fstt;                 /* overall F-statistic */
  double chisq;                /* overall chi-square statistic */
  double lnL;                  /* log-likelihood */
  double ybar, sdy;            /* mean and std. dev. of dependent var. */
  double[ModelSelCriteria.max] criterion;     /* array of model selection statistics */
  double dw, rho;              /* Durbin-Watson stat. and estimated 1st
				    order autocorrelation coefficient */
  ARInfo *arinfo;              /* pointer to struct to hold special info for 
				    autoregressive model */ 
  int errcode;                 /* Error code in case of failure */
  char *name;                  /* for use in GUI */
  char *depvar;                /* name of dependent var in special cases */
  int nparams;                 /* number of named model parameters */
  char **params;               /* for named model parameters */
  version(old) {}
  else {
		long esttime;								 /* it's gint in Gretl; hopefully this is correct */
	}
  int ntests;                  /* number of attached test results */
  ModelTest *tests;            /* attached hypothesis test results */
  Dataset *dataset;            /* for handling models estimated on a sub-sampled portion of the dataset */
  int n_data_items;            /* number of extra data items */
	ModelDataItem **data_items; /* pointer to additional data */
}

Dataset * loadDataset(string filetype="csv")(string filename, bool silent=true) {
	Dataset * result = datainfo_new();
	GretlPrn * p = null;
	if (!silent) {
		p = gretl_print_new(PrintType.stdout, null);
	}
	gretl_read_foreign_data(toUTFz!(char*)(filename), mixin("GretlFileType." ~ filetype), result, p);
	gretl_print_destroy(p);
	return result;
}

DoubleMatrix loadMatrix(string filetype="csv")(string filename, bool silent=true) {
	Dataset * temp = loadDataset!filetype(filename);
	auto result = DoubleMatrix(temp.n, temp.v);
	foreach(col; 0..result.cols) {
		foreach(row; 0..result.rows) {
			result[row, col] = temp.Z[col][row];
		}
	}
	destroy_dataset(temp);
	return result;
}

struct JohansenInfo {}

struct GretlVar {
	int ci;              /* command index (VAR or VECM) */
	int refcount;        /* for saving/deleting */
	int err;             /* error code */
	int neqns;           /* number of equations in system */
	int order;           /* maximum lag order */
	int t1;              /* starting observation */
	int t2;              /* ending observation */
	int T;               /* number of observations */
	int df;              /* T - average coeffs per equation */
	int ifc;             /* equations include a constant (1) or not (0) */
	int ncoeff;          /* total coefficients per equation */
	int * lags;           /* list of specific lags */
	int * ylist;          /* list of stochastic vars */
	int * xlist;          /* list of exogenous variables */
	int * rlist;          /* restricted exogenous variables (VECM only) */
	int detflags;        /* record of automatic deterministic vars added */
	int robust;          /* computing robust std errors? */
	int xcols;           /* full column size of X matrix (VECM special) */
	GretlMatrix * Y;     /* matrix of dependent variables */
	GretlMatrix * X;     /* matrix of independent variables */
	GretlMatrix * B;     /* basic coefficient matrix */
	GretlMatrix * XTX;   /* X'X inverse */
	GretlMatrix * A;     /* augmented coefficient matrix (companion form) */
	GretlMatrix * L;     /* lambda: inverse roots of A(L) polynomial */
	GretlMatrix * E;     /* residuals matrix */
	GretlMatrix * C;     /* augmented Cholesky-decomposed error matrix */
	GretlMatrix * S;     /* cross-equation variance matrix */
	GretlMatrix * F;     /* optional forecast matrix */
	GretlMatrix * ord;   /* optional Cholesky-ordering vector */
	Model** models;      /* pointers to individual equation estimates */
	double * Fvals;       /* hold results of F-tests */
	double * Ivals;       /* hold results of info criteria comparisons */
	double ldet;         /* log-determinant of S */
	double ll;           /* log-likelihood */
	double AIC;          /* Akaike criterion */
	double BIC;          /* Bayesian criterion */
	double HQC;          /* Hannan-Quinn criterion */
	double LR;           /* for likelihood-ratio testing */
	double LB;           /* Ljung-Box (Portmanteau) test statistic */
	int LBs;             /* order for Portmanteau test */
	JohansenInfo * jinfo; /* extra information for VECMs */
	char* name;          /* for use in session management */
};


/*
This function returns the sequence from `val0` to `val1` *inclusive*. It
is inclusive to match R's seq function.
*/

int[] seq(int val0, int val1) {
	int[] result;
	result.reserve(val1-val0+1);
	foreach(val; val0..val1+1) {
		result ~= val;
	}
	return result;
}

extern (C) {
  int gretl_matrix_multiply(const GretlMatrix * a, const GretlMatrix * b, GretlMatrix * c);
  void gretl_matrix_multiply_by_scalar(GretlMatrix * m, double x);
  int gretl_matrix_multiply_mod(const GretlMatrix * a, matmod amod, const GretlMatrix * b, matmod bmod, GretlMatrix * c, matmod cmod);
  int gretl_matrix_cholesky_decomp(GretlMatrix * a);
  int gretl_matrix_kronecker_product(const GretlMatrix * A, const GretlMatrix * B, GretlMatrix * K);
  int gretl_LU_solve(GretlMatrix * a, GretlMatrix * b);
  int gretl_invert_matrix(GretlMatrix * a);
  double gretl_matrix_determinant(GretlMatrix * a, int * err);
  double gretl_matrix_log_determinant(GretlMatrix * a, int * err);
  double gretl_matrix_trace(const GretlMatrix * m);
  void gretl_matrix_raise(GretlMatrix * m, double x);
  void gretl_matrix_free(GretlMatrix * m);
  int gretl_matrix_ols (const GretlMatrix * y, const GretlMatrix * X, GretlMatrix * b, GretlMatrix * vcv, GretlMatrix * uhat, double * s2);
  void gretl_matrix_print(const GretlMatrix * m, const char * msg);
  int gretl_matrix_transpose(GretlMatrix * targ, const GretlMatrix * src);
  GretlMatrix * gretl_matrix_alloc(int rows, int cols);
  double gretl_vector_mean (const GretlMatrix * v);
  double gretl_vector_variance (const GretlMatrix * v);
  
	double gretl_min (int t1, int t2, const double *x);
	double gretl_max (int t1, int t2, const double *x);
	double gretl_sum (int t1, int t2, const double *x);
	double gretl_mean (int t1, int t2, const double *x);
	int gretl_array_quantiles (double *a, int n, double *p, int k);
  double gretl_array_quantile (double *a, int n, double p);
	double gretl_median (int t1, int t2, const double * x);
  double gretl_sst (int t1, int t2, const double *x);
	double gretl_variance (int t1, int t2, const double *x);
	double gretl_stddev (int t1, int t2, const double *x);
	double gretl_covar (int t1, int t2, const double *x, const double *y, int *missing);
	double gretl_corr (int t1, int t2, const double *x, const double *y, int *missing);
	double gretl_skewness (int t1, int t2, const double *x);
	double gretl_kurtosis (int t1, int t2, const double *x);
	int gretl_moments (int t1, int t2, const double *x, const double *wts, double *xbar, double *sd, double *skew, double *kurt, int k);
	
  void gretl_rand_init();
  void gretl_rand_free();
  void gretl_rand_set_seed(uint seed);
  void gretl_rand_set_multi_seed(const GretlMatrix * seed);
  uint gretl_rand_int();
  uint gretl_rand_int_max(uint m);
	int gretl_rand_int_minmax(int *a, int n, int min, int max);
  double gretl_rand_01();
  double gretl_one_snormal();
  void gretl_rand_uniform(double * a, int t1, int t2);
  int gretl_rand_uniform_minmax(double * a, int t1, int t2, double min, double max);
  void gretl_rand_normal(double * a, int t1, int t2);
  int gretl_rand_normal_full(double * a, int t1, int t2, double mean, double sd);
  int gretl_rand_chisq(double * a, int t1, int t2, double v);
  int gretl_rand_student(double * a, int t1, int t2, double v);
  int gretl_rand_F(double * a, int t1, int t2, int v1, int v2);
  int gretl_rand_binomial(double * a, int t1, int t2, int n, double p);
  int gretl_rand_poisson(double * a, int t1, int t2, const double * m, int vec);
  int gretl_rand_weibull(double * a, int t1, int t2, double shape, double scale);
  int gretl_rand_gamma(double * a, int t1, int t2, double shape, double scale);
  double gretl_rand_gamma_one(double shape, double scale);
  int gretl_rand_GED(double * a, int t1, int t2, double nu);
  int gretl_rand_beta(double * x, int t1, int t2, double s1, double s2);
  int gretl_rand_beta_binomial(double * x, int t1, int t2, int n, double s1, double s2);
  GretlMatrix * halton_matrix(int m, int r, int offset, int * err);
  GretlMatrix * inverse_wishart_matrix(const GretlMatrix * S, int v, int * err);
  GretlMatrix * inverse_wishart_sequence(const GretlMatrix * S, int v, int replics, int * err);
  void gretl_rand_set_box_muller(int s);
  int gretl_rand_get_box_muller();
  uint gretl_rand_get_seed();
  
  int * gretl_list_new(int nterms);
  
  Dataset * create_new_dataset(int nvar, int nobs, int markers);
  void destroy_dataset(Dataset * dset);
	int dataset_set_time_series(Dataset * dset, int pd, int yr0, int minor0);
  Dataset * datainfo_new();
  int gretl_read_foreign_data(const char *fname, GretlFileType file_type, Dataset * dset, GretlPrn * prn);
  
  Model * gretl_model_new();
  void gretl_model_free(Model * pmod);
  void clear_model(Model * pmod);
  int printmodel(Model * pmod, const Dataset * dset, int opt, GretlPrn * prn);
	double gretl_model_get_vcv_element(const(Model) * pmod, int i, int j, int np);
	
  Model arma(const(int)* list, const(int)* pqlags, const Dataset * dset, int opt, GretlPrn * prn);
  Model garch(const(int)* list, Dataset * dset, gretlopt opt, GretlPrn * prn);
  Model binary_logit(const int * list, Dataset * dset, gretlopt opt, GretlPrn * prn);
	Model binary_probit(const int *list, Dataset * dset, gretlopt opt, GretlPrn * prn);  
	
	GretlPrn * gretl_print_new(PrintType ptype, int *err);
	GretlPrn * gretl_print_new_with_filename(const char *fname, int *err);
	GretlPrn * gretl_print_new_with_tempfile(int *err);
	GretlPrn * gretl_print_new_with_buffer(char *buf);
  void gretl_print_destroy(GretlPrn * prn);
	
	//~ int adf_test(int order, const int *list, Dataset * dset, 
	      //~ gretlopt opt, GretlPrn *prn);
	GretlVar * gretl_VAR(int order, int * laglist, int * list, 
		      const Dataset *dset, gretlopt opt, GretlPrn * prn, int * err);
	Model * gretl_VAR_get_model(GretlVar * var, int i);
}
