#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>
#include <cmath>
#include "plate.hpp"
#include "cost_function.hpp"


template<typename scalar>
class optimizer {
private:
	plate current_plate;
	cost_function cost;
public:
	optimizer() {
		
	}
	
};