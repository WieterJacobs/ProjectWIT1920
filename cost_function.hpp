#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>
#include <cmath>
#include "plate.hpp"


template<typename scalar>
class cost_function {
private:
	typedef Matrix<scalar, Dynamic, 1> DV;
	typedef Matrix<scalar, 1, Dynamic> DR;
	typedef Array<scalar, Dynamic, Dynamic> DA;
	scalar cost_;
	DR cost_p_;
	void replace_plate(plate<scalar> new_plate) 
	{
		cost_ = cost_f(new_plate.get_T());
		DV lambda;
		//std::cout << cost_T(new_plate.get_T()).transpose() << std::endl;
		//std::cout << new_plate.get_K().transpose() << std::endl;
		lambda = new_plate.solve(new_plate.get_K().transpose(), cost_T(new_plate.get_T()).transpose());
		cost_p_ = -new_plate.get_DK().transpose()*lambda;
		std::cout << cost_p_ << std::endl;
	}
	inline scalar cost_f(DA  T) const 
	{
		return std::pow(T.pow(16).sum(), 1. / 16.);
	}
	inline DV cost_T(DA T) const {
		return T.pow(15) * (std::pow(T.pow(16).sum(), -15. / 16.));
	}
public:
	cost_function(plate<scalar> new_plate)
	{
		replace_plate(new_plate);
	}
	scalar cost() const {
		return cost_;
	}
	DR derivative() {
		return cost_p_;
	}
};