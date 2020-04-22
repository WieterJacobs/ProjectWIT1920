#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "cost_function.hpp"
#include <nlopt.hpp>
#include <iomanip>
#include <chrono>
using namespace Eigen;
int counter = 0;

double g(unsigned n, const double* x, double* grad, void* g_data) {
	double constraint = 0;
	for (unsigned i = 0; i < n; i++) {
		grad[i] = 1. / n;
		constraint += x[i];
	}
	constraint = constraint / n - 0.4;
	std::cout << "distance to constraint: " << -constraint << std::endl;
	return constraint;//zet juist
}
double f(unsigned n, const double* x, double* grad, void* f_data) {	
	cost_function<double>* cost = static_cast<cost_function<double>*>(f_data);
	cost->update_cost(x);
	for (int i = 0; i < n; i++) {
		grad[i] = cost->get_cost_v()[i];
	}
	double f_val= cost->get_cost();
	counter++;
	std::cout <<counter << "	| " <<"Objective function value: " << f_val << "	|	";
	return(f_val);//cost->get_cost());
}
int main() {

	std::cout << std::setprecision(6);
	int n = 20;
	double size = 0.01;
	
	ArrayXXd v;
	v = ArrayXXd  ::Constant(n, n, 0.35);
	double* v_data = v.data();
	double p = 1;
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_MMA, n*n);
	nlopt_set_lower_bounds1(opt, 0.);
	nlopt_set_upper_bounds1(opt, 1.);
	cost_function<double>* cost = new cost_function<double>(n, size, p);
	nlopt_set_min_objective(opt, f, cost);
	nlopt_add_inequality_constraint(opt, g, NULL, 0);
	nlopt_set_xtol_rel(opt, 1e-2);
	nlopt_set_ftol_rel(opt, 1e-4);
	double minf;

	auto start = std::chrono::high_resolution_clock::now();

	std::cout << "---------" << "p = " << p << "--------" << std::endl;
	nlopt_optimize(opt, v_data, &minf);

	p = 2; 
	std::cout << "---------" << "p = " << p << "--------" << std::endl;
	cost->set_p(p);
	nlopt_set_min_objective(opt, f, cost);
	nlopt_set_xtol_rel(opt, 1e-4);
	nlopt_set_ftol_rel(opt, 1e-7);
	nlopt_optimize(opt, v_data, &minf);
	
	p = 3;
	std::cout << "---------" << "p = " << p << "--------" << std::endl;
	cost->set_p(p);
	nlopt_set_min_objective(opt, f, cost);
	nlopt_set_xtol_rel(opt, 1e-5);
	nlopt_set_ftol_rel(opt, 1e-8);
	nlopt_optimize(opt, v_data, &minf);
	auto stop = std::chrono::high_resolution_clock::now();
	nlopt_destroy(opt);
	delete cost;
	
	ArrayXXd v_tl = Map<ArrayXXd>(v_data, n, n);
	ArrayXXd v_bl = v_tl.block(0, 0, n - 1, n).colwise().reverse();
	ArrayXXd v_tr = v_tl.block(0, 0, n, n - 1).rowwise().reverse();
	ArrayXXd v_br = v_tl.block(0, 0, n-1, n - 1).reverse();
	ArrayXXd result(2*n - 1, 2*n - 1);
	delete[] v_data;
	//delete v_data;
	result.block(0,0,n,n)<<v_tl;
	result.block(0, n, n, n - 1) << v_tr;
	result.block(n, 0, n - 1, n) << v_bl;
	result.block(n, n, n - 1, n - 1) << v_br;
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	std::cout << "minimum found in " << duration.count()*1e-6 << " seconds" << std::endl;
	
	//test stuff
	n = 20;
	ArrayXXd v_test = ArrayXXd::Constant(n, n, 0.4);
	v_test.setRandom();
	v_test = (v_test + 1) / 2;
	double* v_test_data = v_test.data();
	p = 3;
	size = 0.01;
	cost_function<double> cost_f = cost_function<double>(n, size, p);
	cost_f.update_cost(v_test_data);
	ArrayXXd cost_v = Map<ArrayXXd>(cost_f.get_cost_v(), n, n);
	plate<double> heated_plate = plate<double>(v_test,size,p);
	MatrixXd K = heated_plate.get_K();
	MatrixXd DK = heated_plate.get_DK();
	MatrixXd T = heated_plate.get_T();
	delete[] v_test_data;

	//simulation speed test
	n = 300;
	ArrayXXd v_test_2 = ArrayXXd::Constant(n, n, 0.4);
	v_test_2.setRandom();
	v_test_2 = (v_test_2 + 1) / 2;
	plate<double> heated_plate_2 = plate<double>(v_test_2, size, p);
	duration = heated_plate_2.simulationTime();
	std::cout << "Simulation speed test of size " << n << " completed in " << duration.count()*1e-6 << " seconds." << std::endl;

	//writing to file
	std::ofstream myfile;
	myfile.precision(17);
	myfile.open("result.txt");
	myfile << result << std::endl;
	myfile.close();
	myfile.open("v_matrix.txt");
	myfile << v_test << std::endl;
	myfile.close();
	myfile.open("K_matrix.txt");
	myfile << MatrixXd(K) << std::endl;
	myfile.close();
	myfile.open("DK_matrix.txt");
	myfile << MatrixXd(DK) << std::endl;
	myfile.close();
	myfile.open("T_matrix.txt");
	myfile << T << std::endl;
	myfile.close();
	myfile.open("cost_v_matrix.txt");
	myfile << cost_v << std::endl;
	myfile.close();
	std::cout << "end" << std::endl;
	return 0;
}