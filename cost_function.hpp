#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>
#include <cmath>
#include "plate.hpp"
#include <nlopt.hpp>





template<typename scalar>
class cost_function {
private:
	typedef Matrix<scalar, Dynamic, 1> DV;
	typedef Matrix<scalar, Dynamic, Dynamic> DM;
	typedef Matrix<scalar, 1, Dynamic> DR;
	typedef Array<scalar, Dynamic, Dynamic> DA;

	int n_;
	scalar size_;
	scalar p_;

	scalar cost_;
	DM cost_v_;
	void adjoint(plate<scalar> new_plate)
	{
		cost_ = cost_f(new_plate.get_T());
		DV lambda;
		lambda = new_plate.solve(new_plate.get_K().transpose(), cost_T(new_plate.get_T()).transpose());
		cost_v_ = -new_plate.get_DK().transpose()*lambda;
		DA  c_v = new_plate.conductivity_v();
		cost_v_.resizeLike(c_v);
		cost_v_ = cost_v_.array() * c_v;
	}
	inline scalar cost_f(DA  T) const 
	{
		return T.sum();
		//return std::pow(T.pow(16).sum(), 1. / 16.);
	}
	inline DV cost_T(DA T) const {
		return T.pow(0);
		//return T.pow(15) * (std::pow(T.pow(16).sum(), -15. / 16.));
	}
	
public:
	cost_function(int n, scalar size, scalar p) :
		n_(n),
		size_(size),
		p_(p) {}
	void update_cost(const double* v_data) {
		const DA v = Map<const DA>(v_data, n_, n_);
		plate<scalar> new_plate = plate<scalar>(v, size_, p_);
		adjoint(new_plate);
	}
	double get_cost() { return cost_; }
	double* get_cost_v() {return cost_v_.data(); }
};