#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "cost_function.hpp"


int main() {
	ArrayXXd v;
	v = ArrayXXd  ::Constant(5, 5, 0.4);
	plate<double> heated_plate = plate<double>(v, 0.01);

	cost_function<double> cost_f = cost_function<double>(heated_plate);
	//std::cout << cost_f.cost() << std::endl;
	//std::cout << cost_f.derivative() << std::endl;
	MatrixXd K = heated_plate.get_K();
	MatrixXd DK = heated_plate.get_DK();
	MatrixXd T = heated_plate.get_T();
	std::ofstream myfile;
	myfile.open("K_matrix.txt");
	myfile << K << std::endl;
	myfile.close();
	myfile.open("DK_matrix.txt");
	myfile << DK << std::endl;
	myfile.close();
	myfile.open("T_matrix.txt");
	myfile << T << std::endl;
	myfile.close();
	return 0;
}