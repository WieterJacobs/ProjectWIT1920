#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>
#include <cmath>
#include <chrono>

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
		int ldir_ =round(.6 * (n_-2));
		int udir_ = n_-1;
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
				if (i >= ldir_ && i<udir_)//vraag over stellen aan wieter
				{
					coefficients.push_back(T(pos(i, 0), pos(i, 0), 1));
					Q_(pos(i, 0)) = 293;
					coefficients.push_back(T(pos(i, n_ - 1), pos(i, n_ - 1), 1));
					//Q_(pos(i, n_-1)) = 293;
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
			Matrix<scalar,4,4> dK;
			Matrix<scalar,4,1> temp;
			scalar north;
			scalar east;
			scalar south;
			scalar west;
			for (int i = 0; i < n_ - 1; ++i)
			{
				for (int j = 0; j < n_ - 1; ++j)
				{
					if (i>0) {north = dmeanh(k(i, j), k(i - 1, j)); }
					else { north = 0; }
					if (j<n_-2) { east = dmeanh(k(i, j), k(i, j + 1));}
					else { east = 0; }
					if (i<n_-2) { south = dmeanh(k(i, j), k(i + 1, j)); }
					else { south = 0; }
					if (j>0) { west = dmeanh(k(i, j), k(i, j - 1)); } 
					else { west = 0; } 
					dK << south+west, -south, 0, -west,
						-south, east+south, -east, 0,
						0, -east, east+north, -north,
						-west, 0, -north, west+north;
					temp << T_(pos(i + 1, j)), T_(pos(i + 1, j + 1)), T_(pos(i, j + 1)), T_(pos(i, j));
					temp = dK * temp;

					if (j != 0) {
						if (i != 0) {
							coefficients.push_back(T(pos(i, j), posv(i, j), temp(3)));
						}
						if (i != n_-2 ) {
							coefficients.push_back(T(pos(i + 1, j), posv(i, j), temp(0)));
						}
					}
					if (j != n_-2) {
						if (i != 0) {
							coefficients.push_back(T(pos(i, j + 1), posv(i, j), temp(2)));
						}
						if (i != n_-2) {
							coefficients.push_back(T(pos(i + 1, j + 1), posv(i, j), temp(1)));
						}
					}
					/*
					if (!(j == 0 && i > ldir_ -2)) {
						coefficients.push_back(T(pos(i + 1, j), posv(i, j), temp(0)));
					}
					coefficients.push_back(T(pos(i + 1, j + 1), posv(i, j), temp(1)));
					coefficients.push_back(T(pos(i, j + 1), posv(i, j), temp(2)));
					if (!(j == 0 && i >= ldir_ -2)) {
						coefficients.push_back(T(pos(i, j), posv(i, j), temp(3)));
					}
					*/
				}
			}
			DK_.setFromTriplets(coefficients.begin(), coefficients.end());
		}
	public:
		inline DM conductivity_v() const &
		{

			return p_ * k_metal_ * v_.pow(p_ - 1) - p_ * k_plastic_ * v_.pow(p_ - 1);

		}

		plate(const DA& v, scalar const size, scalar p) :
			size_(size),
			deltaz_(.001),
			n_(v.rows()+1),
			deltax_(size_ / n_),
			N_(n_* n_),
			p_(p),
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
		auto simulationTime() {
			auto start = std::chrono::high_resolution_clock::now();
			set_K();
			T_ = solve(K_, Q_);
			auto stop = std::chrono::high_resolution_clock::now();
			return std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		}
};