import std.conv, std.exception, std.math, std.stdio;
import gretl.base, gretl.matrix, gretl.random;

int adf_test(int order, int[] list, Dataset * dset, 
	      gretlopt opt, GretlPrn * prn) {
  int save_t1 = dset.t1;
  int save_t2 = dset.t2;
    //~ int panelmode;
  int err;

    /* GLS incompatible with no const, quadratic trend or seasonals */
    // Figure this out, for now ignore the options and go with the usual case
    //~ err = incompatible_options(opt, OPT_G | OPT_N | OPT_R);
    //~ if (!err) {
			//~ err = incompatible_options(opt, OPT_D | OPT_G);
    //~ }

    //~ if (!err && (opt & OPT_G)) {
		/* under GLS, have to choose between cases */
			//~ err = incompatible_options(opt, OPT_C | OPT_T);
    //~ }

    //~ panelmode = multi_unit_panel_sample(dset);

    //~ if (panelmode) {
			//~ err = panel_DF_test(list[1], order, dset, opt, prn);
    //~ } else {
			/* regular time series case */
			int i, v;
			int[2] vlist = [1, 0];
			auto ainfo = AdfInfo(0);
			ainfo.niv = 1;

			foreach(ii; 0..list.length-1) {
				ainfo.v = vlist[1] = list[ii+1];
				ainfo.order = order;
				vlist[1] = v = list[ii+1];
				err = list_adjust_sample(vlist.ptr, &(dset.t1), &(dset.t2), dset, null);
				
				// Set order to -1 to use a default value
				if (!err && order == -1) {
					/* default to L_{12}: see G. W. Schwert, "Tests for Unit Roots:
						 A Monte Carlo Investigation", Journal of Business and
						 Economic Statistics, 7(2), 1989, pp. 5-17. Note that at
						 some points Ng uses floor(T/100.0) in the following
						 expression, which can give a lower max order.
					*/
					int T = dset.t2 - dset.t1 + 1;
					double tmp = 12.0 * pow(T/100.0, 0.25);
					ainfo.order = tmp.to!int;
				}
				writeln("err ", err);
				if (!err) {
					writeln("Before ", err);
					err = real_adf_test(&ainfo, dset, gretlopt.none, null);
					writeln("After ", err);
				}
				dset.t1 = save_t1;
				dset.t2 = save_t2;
			}
    //~ }
    dset.t1 = save_t1;
    dset.t2 = save_t2;
    return err;
}

void main() {
	randInit();
	setSeed(200);
	Dataset * ds = create_new_dataset(1, 100, 0);
	scope(exit) { destroy_dataset(ds); }
	foreach(ii; 0..100) {
		(*ds)[ii,0] = rnorm();
	}
	dataset_set_time_series(ds, 1, 1, 1);
	writeln(adf_test(2, [1, 0], ds, gretlopt.none, null));
}

struct AdfInfo {
	int v;           /* ID number of series to test (in/out) */
	int order;       /* lag order for ADF (in/out) */
	int kmax;        /* max. order (for testing down) */
	int altv;        /* ID of modified series (detrended) */
	int niv;         /* number of (co-)integrated vars (Engle-Granger) */
	AdfFlags flags;  /* bitflags: see above */
	DetCode det;     /* code for deterministics */
	int nseas;       /* number of seasonal dummies */
	int T;           /* number of obs used in test */
	int df;          /* degrees of freedom, test regression */
	double b0;       /* coefficient on lagged level */
	double tau;      /* test statistic */
	double pval;     /* p-value of test stat */
	int * list;       /* regression list */
	int * slist;      /* list of seasonal dummies, if applicable */
	char *vname; /* name of series tested */
	GretlMatrix * g; /* GLS coefficients (if applicable) */
};

enum AdfFlags {
	egtest   = 1 << 0, /* doing Engle-Granger test */
	egresids = 1 << 1, /* final stage of the above */
	panel     = 1 << 2, /* working on panel data */
	olsfirst = 1 << 3  /* Perron-Qu, 2007 */
}

enum DetCode {
	noconst = 1,
	urconst,
	trend,
	quad,
	max
} 

int real_adf_test(AdfInfo * ainfo, Dataset * dset,
			  gretlopt opt, GretlPrn * prn) {
    Model dfmod;
    gretlopt eg_opt = gretlopt.none;
    gretlopt df_mod_opt = (gretlopt.a | gretlopt.z);
    int * biglist = null;
    int orig_nvars = dset.v;
    int blurb_done = 0;
    int test_down = 0;
    int test_num = 0;
    int i, err;

    /* (most of) this may have been done already
       but it won't hurt to check here */
    //~ err = check_adf_options(opt);
    //~ if (err) {
	//~ return err;
    //~ }

    //~ /* safety-first initializations */
    ainfo.nseas = ainfo.kmax = ainfo.altv = 0;
    ainfo.vname = dset.varname[ainfo.v];
    ainfo.list = ainfo.slist = null;
    ainfo.det = 0;

	return 0;
}

extern(C) {
	int list_adjust_sample(const int * list, int * t1, int * t2, 
			const Dataset * dset, int * nmiss);
}
