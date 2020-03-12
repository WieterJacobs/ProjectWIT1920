#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
//#include "cost_function.hpp"
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"


using namespace alglib;
void function1_grad(const real_1d_array& x, double& func, real_1d_array& grad, void* ptr)
{
    //
    // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
    // and its derivatives df/d0 and df/dx1
    //
    func = 100 * pow(x[0] + 3, 4) + pow(x[1] - 3, 4);
    grad[0] = 400 * pow(x[0] + 3, 3);
    grad[1] = 4 * pow(x[1] - 3, 3);
}

int main(int argc, char** argv)
{
    //
    // This example demonstrates minimization of
    //
    //     f(x,y) = 100*(x+3)^4+(y-3)^4
    //
    // subject to inequality constraints
    //
    // * x>=2 (posed as general linear constraint),
    // * x+y>=6
    //
    // using BLEIC optimizer with
    // * initial point x=[0,0]
    // * unit scale being set for all variables (see minbleicsetscale for more info)
    // * stopping criteria set to "terminate after short enough step"
    // * OptGuard integrity check being used to check problem statement
    //   for some common errors like nonsmoothness or bad analytic gradient
    //
    // First, we create optimizer object and tune its properties:
    // * set linear constraints
    // * set variable scales
    // * set stopping criteria
    //
    real_1d_array x = "[5,5]";
    real_1d_array s = "[1,1]";
    real_2d_array c = "[[1,0,2],[1,1,6]]";
    integer_1d_array ct = "[1,1]";
    minbleicstate state;
    double epsg = 0;
    double epsf = 0;
    double epsx = 0.000001;
    ae_int_t maxits = 0;

    minbleiccreate(x, state);
    minbleicsetlc(state, c, ct);
    minbleicsetscale(state, s);
    minbleicsetcond(state, epsg, epsf, epsx, maxits);

    //
    // Then we activate OptGuard integrity checking.
    //
    // OptGuard monitor helps to catch common coding and problem statement
    // issues, like:
    // * discontinuity of the target function (C0 continuity violation)
    // * nonsmoothness of the target function (C1 continuity violation)
    // * erroneous analytic gradient, i.e. one inconsistent with actual
    //   change in the target/constraints
    //
    // OptGuard is essential for early prototyping stages because such
    // problems often result in premature termination of the optimizer
    // which is really hard to distinguish from the correct termination.
    //
    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
    //            DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
    //
    //            Other OptGuard checks add moderate overhead, but anyway
    //            it is better to turn them off when they are not needed.
    //
    minbleicoptguardsmoothness(state);
    minbleicoptguardgradient(state, 0.001);

    //
    // Optimize and evaluate results
    //
    minbleicreport rep;
    alglib::minbleicoptimize(state, function1_grad);
    minbleicresults(state, x, rep);
    printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
    printf("%s\n", x.tostring(2).c_str()); // EXPECTED: [2,4]

    //
    // Check that OptGuard did not report errors
    //
    // NOTE: want to test OptGuard? Try breaking the gradient - say, add
    //       1.0 to some of its components.
    //
    optguardreport ogrep;
    minbleicoptguardresults(state, ogrep);
    printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
    printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
    printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
    return 0;
}


/*
int main() {
	
	ArrayXXd v;
	v = ArrayXXd  ::Constant(5, 5, 0.4);
	plate<double> heated_plate = plate<double>(v, 0.01,1);

	//cost_function<double> cost_f = cost_function<double>(v,0.01,1);
	//std::cout << cost_f.cost() << std::endl;
	//std::cout << cost_f.derivative() << std::endl;
	MatrixXd K = heated_plate.get_K();
	MatrixXd DK = heated_plate.get_DK();
	MatrixXd T = heated_plate.get_T();
	//MatrixXd cost_p = cost_f.derivative();
	//std::cout.precision(17);
	//std::cout << T << std::endl;
	std::ofstream myfile;
	
	myfile.open("K_matrix.txt");
	myfile.precision(17);
	myfile << K << std::endl;
	myfile.close();
	myfile.open("DK_matrix.txt");
	myfile << DK << std::endl;
	myfile.close();
	myfile.open("T_matrix.txt");
	myfile << T << std::endl;
	myfile.close();
	myfile.open("cost_p_matrix.txt");
	//myfile << cost_p << std::endl;
	myfile.close();
	return 0;
	

}*/