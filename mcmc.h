#ifndef _mcmc_h_
#define _mcmc_h_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include <pthread.h>
#include <csignal>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#define rho_crit 6.361025188e-8
#define PI M_PI

//Properties of the C-M relation for halo prior
#define zero 0.83
#define slope -0.098
#define hubble 0.72
#define sigc 0.105

//Size of scale factor decline constant
#define scale_constant 10000.0

using namespace std;

void           CopyExtraParameters(vector<double> &Extra_Params);
vector<string> FindParameters(string &line_of_data);
double         GenGalaxyAndGetChiSquare();
double         GetConcentration(const double mass);
void           signal_handler(int signum);
void           WriteDataFiles(vector<double> &param);

#endif
