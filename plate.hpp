#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>

using namespace Eigen;

template<typename scalar>
	class plate {
		
	private:
		typedef Triplet<scalar> T;
		typedef std::vector<scalar> V;
		typedef std::vector<T> TV;
		typedef SparseMatrix<scalar> SM;
		typedef Array<scalar, Dynamic, Dynamic> DA;
		typedef Matrix<scalar, Dynamic, 1> RHS;

		scalar k_metal_ = 65.;
		scalar k_plastic_ = .2;
		scalar size_;
		int n_;
		scalar deltax_;
		int N_;
		SM F_;
		RHS Q_;
		scalar ldir_ = floor(.3 * n_);
		scalar udir_ = ceil(.7 * n_);
		inline int pos(int i, int j)
		{
			return(i + n_ * j);
		}
	public:
		plate(DA& v, scalar const size) :
			size_(size),
			n_(v.rows()),
			deltax_(size_ / n_),
			N_(n_* n_),
			F_(N_, N_)
		{
			DA k(n_,n_);
			Q_ = RHS::Constant(N_,0);
			k = k_metal_ * v + k_plastic_ * (1 - v);
			TV coefficients;
			V w(4);
			for (int i = 1; i < n_ - 1; ++i)
			{
				for (int j = 1; j < n_ - 1; ++j)
				{
					Q_(pos(i, j)) = 2. * deltax_*deltax_*10e5;
					w[0] = (k(i - 1, j) + k(i, j)) / 2;
					w[1] = (k(i + 1, j) + k(i, j)) / 2;
					w[2] = (k(i, j - 1) + k(i, j)) / 2;
					w[3] = (k(i, j + 1) + k(i, j)) / 2;

					coefficients.push_back(T(pos(i, j), pos(i - 1, j), -w[0]));
					coefficients.push_back(T(pos(i, j), pos(i + 1, j), -w[1]));
					coefficients.push_back(T(pos(i, j), pos(i, j - 1), -w[2]));
					coefficients.push_back(T(pos(i, j), pos(i, j + 1), -w[3]));
					coefficients.push_back(T(pos(i, j), pos(i, j), std::accumulate(w.begin(), w.end(), 0.0)));
				}
			}
			for (int j = 0; j < n_; ++j)
			{
				coefficients.push_back(T(pos(0, j), pos(0, j), 1));
				coefficients.push_back(T(pos(0, j), pos(1, j), -1));
				coefficients.push_back(T(pos(n_ - 1, j), pos(n_ - 1, j), 1));
				coefficients.push_back(T(pos(n_ - 1, j), pos(n_ - 2, j), -1));
			}
			for (int i = 1; i < n_ - 1; ++i)
			{
				if (i > ldir_&& i < udir_)
				{
					coefficients.push_back(T(pos(i, 0), pos(i, 0), 1));
					Q_(pos(i, 0)) = 293;
					coefficients.push_back(T(pos(i, n_ - 1), pos(i, n_ - 1), 1));
					Q_(pos(i, n_ - 1)) = 293;
				}
				else
				{
					coefficients.push_back(T(pos(i, 0), pos(i, 0), 1));
					coefficients.push_back(T(pos(i, 0), pos(i, 1), -1));
					coefficients.push_back(T(pos(i, n_ - 1), pos(i, n_ - 1), 1));
					coefficients.push_back(T(pos(i, n_ - 1), pos(i, n_ - 2), -1));
				}
			}
			F_.setFromTriplets(coefficients.begin(), coefficients.end());
			//std::cout << MatrixXd(F_) << std::endl;
		}
		DA temperature() const 
		{
			SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;//checken of dit wel een goede keuze van solver is
			solver.analyzePattern(F_);
			solver.factorize(F_);
			RHS x(N_);
			DA y(N_, 1);
			x= solver.solve(Q_);
			y = x;//lijkt niet nodige copy maar ik kreeg het niet op een andere manier gedaan.
			y.resize(n_, n_);
			return(y);
		}
};