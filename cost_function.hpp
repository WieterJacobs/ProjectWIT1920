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
	plate<scalar> current_plate_;
	void adjoint()
	{
		cost_ = cost_f(current_plate_.get_T());
		DV lambda;
		lambda = current_plate_.solve(current_plate_.get_K().transpose(), cost_T(current_plate_.get_T()).transpose());
		cost_p_ = -current_plate_.get_DK().transpose()*lambda;
		DA  c_v = current_plate_.conductivity_v();
		c_v.resizeLike(cost_p_);
		cost_p_ = cost_p_.array() * c_v;
	}
	inline scalar cost_f(DA  T) const 
	{
		return std::pow(T.pow(16).sum(), 1. / 16.);
	}
	inline DV cost_T(DA T) const {
		return T.pow(15) * (std::pow(T.pow(16).sum(), -15. / 16.));
	}
public:

	cost_function(double size, double p)
	{
	}
	double EvaluateWithGradient(const ) const {

	}
};