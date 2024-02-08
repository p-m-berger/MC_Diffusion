
DIFF_EXTERN void calc_diff(double *, double, double, double);
DIFF_EXTERN int calc_dt(double, double *);

DIFF_EXTERN int load_input(char *, double *, double *, double *);
DIFF_EXTERN void write_step(double, int);
DIFF_EXTERN void write_delts(double);
DIFF_EXTERN void write_mins(char *, int);
DIFF_EXTERN void out_d(double **, int, double);
DIFF_EXTERN void panic(char *);

DIFF_EXTERN void run_model(int, node *, double *, double *, double *, bool *, bool);
DIFF_EXTERN void initc(int, node *, int, double, bool, double *, double *);