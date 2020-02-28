#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>
#include <cmath>

using namespace Eigen;

template<typename scalar>
	class plate {
		
	private:
		typedef Triplet<scalar> T;
		typedef std::vector<scalar> V;
		typedef std::vector<T> TV;
		typedef SparseMatrix<scalar> SM;
		typedef Array<scalar, Dynamic, Dynamic> DA;
		typedef Matrix<scalar, Dynamic, Dynamic> DM;
		typedef Matrix<scalar, Dynamic, 1> DV;

		scalar k_metal_ = 65.;
		scalar k_plastic_ = 0.2;
		scalar size_;
		scalar deltaz_;
		int n_;
		scalar deltax_;
		int N_;
		mutable scalar p_;
		DA v_;
		scalar ldir_ = (.6 * n_);
		mutable DV T_;
		DV Q_;
		SM K_;
		SM DK_;
		inline int pos(int i, int j) const
		{
			return(i + n_ * j);
		}
		inline int posv(int i, int j) const
		{
			return(i + (n_ - 1) * j);
		}
		inline scalar meanh(scalar k1, scalar k2) const {
			return(2 * k1 * k2 / (k1 + k2));
		}
		inline scalar dmeanh(scalar k1, scalar k2) const {
			return(2 * k2 * k2 / ((k1 + k2) * (k1 + k2)));
		}
		inline DA conductivity(DA v) const
		{
			return( k_metal_ * v.pow(p_) + k_plastic_ * (1 - v.pow(p_)));
		}
		inline scalar conductivity_v() const
		{
			//todo
		}
		
		void set_K() {
			DA k = conductivity(v_);
			Q_ = DV::Constant(N_, 0);
			TV coefficients;
			V w(4);
			for (int i = 1; i < n_ - 1; ++i)
			{
				for (int j = 1; j < n_ - 1; ++j)
				{
					Q_(pos(i, j)) = 2. * deltax_ * deltax_ / (4 * size_ * size_ * deltaz_);
					w[0] = meanh(k(i - 1, j - 1), k(i - 1, j));
					w[1] = meanh(k(i, j - 1), k(i, j));
					w[2] = meanh(k(i - 1, j - 1), k(i, j - 1));
					w[3] = meanh(k(i - 1, j), k(i, j));

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
				if (i >= ldir_)
				{
					coefficients.push_back(T(pos(i, 0), pos(i, 0), 1));
					Q_(pos(i, 0)) = 293;
					coefficients.push_back(T(pos(i, n_ - 1), pos(i, n_ - 1), 1));
					coefficients.push_back(T(pos(i, n_ - 1), pos(i, n_ - 2), -1));
				}
				else
				{
					coefficients.push_back(T(pos(i, 0), pos(i, 0), 1));
					coefficients.push_back(T(pos(i, 0), pos(i, 1), -1));
					coefficients.push_back(T(pos(i, n_ - 1), pos(i, n_ - 1), 1));
					coefficients.push_back(T(pos(i, n_ - 1), pos(i, n_ - 2), -1));
				}
			}
			K_.setFromTriplets(coefficients.begin(), coefficients.end());
		}
		void set_DK()
		{
			DA k = conductivity(v_);
			TV coefficients;
			Matrix4d dK;//Verander hier misschien Matrix4d voor elegantie (geen gebruik scalar)
			Vector4d temp;//same
			for (int i = 1; i < n_ - 2; ++i)
			{
				for (int j = 1; j < n_ - 2; ++j)
				{
					dK << dmeanh(k(i, j), k(i + 1, j)) + dmeanh(k(i, j), k(i, j - 1)), -dmeanh(k(i, j), k(i + 1, j)), 0, -dmeanh(k(i, j), k(i, j - 1)),
						-dmeanh(k(i, j), k(i + 1, j)), dmeanh(k(i, j), k(i + 1, j)) + dmeanh(k(i, j), k(i, j + 1)), -dmeanh(k(i, j), k(i, j + 1)), 0,
						0, -dmeanh(k(i, j), k(i, j + 1)), dmeanh(k(i, j), k(i, j + 1)) + dmeanh(k(i, j), k(i - 1, j)), -dmeanh(k(i, j), k(i - 1, j)),
						-dmeanh(k(i, j), k(i, j - 1)), 0, -dmeanh(k(i, j), k(i - 1, j)), dmeanh(k(i, j), k(i, j - 1)) + dmeanh(k(i, j), k(i - 1, j));
					temp << T_(pos(i + 1, j)), T_(pos(i + 1, j + 1)), T_(pos(i, j + 1)), T_(pos(i, j));
					temp = dK * temp;
					coefficients.push_back(T(pos(i + 1, j), posv(i, j), temp(0)));
					coefficients.push_back(T(pos(i + 1, j + 1), posv(i, j), temp(1)));
					coefficients.push_back(T(pos(i, j + 1), posv(i, j), temp(2)));
					coefficients.push_back(T(pos(i, j), posv(i, j), temp(3)));
				}
			}
			DK_.setFromTriplets(coefficients.begin(), coefficients.end());
		}
	public:
		plate(DA& v, scalar const size) :
			size_(size),
			deltaz_(.001),
			n_(v.rows()),
			deltax_(size_ / n_),
			N_(n_* n_),
			p_(1),
			v_(v),
			Q_(N_),
			K_(N_,N_),
			DK_(N_, (n_ - 1)* (n_ - 1))
		{
			set_K();
			T_ = solve(K_, Q_);
			set_DK();
		}
		SM get_K() 
		{
			return K_;
		}
		SM get_DK()
		{
			return DK_;
		}
		DA get_T() const &{
			return(T_);
		}
		void set_p(scalar p) const {
			p_ = p;
		}
		scalar get_p() const {
			return(p_);
		}
		scalar get_n() const {
			return n_;
		}
		scalar get_N() const {
			return N_;
		}
		DV solve(SM lhs, DV rhs) const
		{
			SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;//checken of dit wel een goede keuze van solver is
			solver.analyzePattern(lhs);
			solver.factorize(lhs);
			DV x(rhs.rows());
			x = solver.solve(rhs);
			return(x);
		}
};