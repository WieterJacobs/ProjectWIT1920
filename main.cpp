#include <Eigen/Dense>
#include "flux_matrix.hpp"

int main() {
	ArrayXXd v;
	v = ArrayXXd::Constant(5, 5, 1);
	flux_matrix<double> flux = flux_matrix<double>(v, 0.1);
}
