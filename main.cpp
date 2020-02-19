#include <Eigen/Dense>
#include <iostream>
#include "plate.hpp"


int main() {
	ArrayXXd v;
	ArrayXXd T;
	v = ArrayXXd  ::Constant(20, 20, 1);
	plate<double> heated_plate = plate<double>(v, 0.1);
	T = heated_plate.temperature();
	std::cout << T << std::endl;
}