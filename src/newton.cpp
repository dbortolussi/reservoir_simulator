#pragma once
#include <vector>
#include <cmath>
#include <Eigen/Sparse>
#include "model.cpp"
#include <iostream>

struct NewtonReport {
	bool converged;
	int iters;
};

inline Eigen::VectorXd assembleResidual(const Model& m, const Eigen::VectorXd& x,
	const std::vector<double>& Sw_old, double dt) {
	
	const Grid& g = m.grid;
	const Fluid& fl = m.fluid;
	int N = g.numCells();
	double V = g.cellVolume();
	Eigen::VectorXd R = Eigen::VectorXd::Zero(2 * N);

	for (int c = 0;c < N;++c) {
		double pc = x(2 * c);
		double Sc = x(2 * c + 1);
		const Well* w = m.wellAt(c);
		if (w) {
			R(2 * c) = pc - w->p_bhp;
			if (w->type == Well::Injector) {
				R(2 * c + 1) = Sc - (1.0 - fl.kr.Sor);
			}
			else
				R(2 * c + 1) = Sc - Sw_old[c];
			continue;
		}
		R(2 * c + 1) += (m.phi * V / dt) * (Sc - Sw_old[c]);
	}

	for (int c = 0; c < N;++c) {
		for (const Neighbor& f : g.neighbors(c)) {
			int j = f.id;
			if (j < c) continue;
			double T = transmissibility(f, m.k(c), m.k(j));
			double dp = x(2*c) - x(2*j);
			int up = (dp >= 0.0) ? c : j;
			double Sup = x(2 * up + 1);
			double F = T * fl.lambda_t(Sup) * dp;
			double Fw = fl.fw(Sup) * F;
			
			const Well* wc = m.wellAt(c);
			const Well* wj = m.wellAt(j);
			if (!wc) R(2 * c) += F;
			if (!wj) R(2 * j) -= F;
			if (wc && wc->type == Well::Injector) {
				if (!wj || wj->type != Well::Injector) R(2 * j + 1) -= Fw;
			}
			else if (wj && wj->type == Well::Injector) {
				if (!wc || wc->type != Well::Injector) R(2 * c + 1) += Fw;
			}
			else {
				if (!wc) R(2 * c + 1) += Fw;
				if (!wj) R(2 * j + 1) -= Fw;
			}
		}
	}
	return R;
}

inline Eigen::SparseMatrix<double> assembleJacobian(const Model& m, const Eigen::VectorXd& x, double dt) {

	const Grid& g = m.grid;
	const Fluid& fl = m.fluid;
	int N = g.numCells();
	double V = g.cellVolume();
	std::vector<Eigen::Triplet<double>> trip;
	trip.reserve(N * 12);

	for (int c = 0; c < N; ++c) {
		const Well* w = m.wellAt(c);
		if (w) {
			trip.emplace_back(2 * c, 2 * c, 1.0);
			trip.emplace_back(2 * c + 1, 2 * c + 1, 1.0);
			continue;
		}
		trip.emplace_back(2 * c + 1, 2 * c + 1, m.phi * V / dt);
	}

	for (int c = 0;c < N;++c) {
		for (const Neighbor& f : g.neighbors(c)) {
			int j = f.id;
			if (j < c) continue;
			double dp = x(2 * c) - x(2 * j);
			int up = (dp >= 0.0) ? c : j;
			int upS = 2 * up + 1;
			double Sup = x(upS);
			double T = transmissibility(f, m.k(c), m.k(j));
			double lam = fl.lambda_t(Sup);
			double dlam = fl.dlambda_t_dSw(Sup);
			double fw = fl.fw(Sup);
			double dfw = fl.dfw_dSw(Sup);

			double dR_dpi = T * lam;
			double dR_dpj = -T * lam;
			double dR_dS = T * dlam * dp;
			double dRw_dpi = fw * dR_dpi;
			double dRw_dpj = fw * dR_dpj;
			double dRw_dS = T * dp * (dfw * lam + fw * dlam);

			const Well* wc = m.wellAt(c);
			const Well* wj = m.wellAt(j);
			bool skipP_c = (wc != nullptr);
			bool skipP_j = (wj != nullptr);
			bool skipS_c = (wc && wc->type == Well::Injector);
			bool skipS_j = (wj && wj->type == Well::Injector);

			if (!skipP_c) {
				trip.emplace_back(2 * c, 2 * c, dR_dpi);
				trip.emplace_back(2 * c, 2 * j, dR_dpj);
				trip.emplace_back(2 * c, upS, dR_dS);
			}
			if (!skipP_j) {
				trip.emplace_back(2 * j, 2 * c, -dR_dpi);
				trip.emplace_back(2 * j, 2 * j, -dR_dpj);
				trip.emplace_back(2 * j, upS, -dR_dS);
			}

			if (wc && wc->type == Well::Injector) {
				if (!wj || wj->type != Well::Injector) {
					trip.emplace_back(2 * j + 1, 2 * c, -dRw_dpi);
					trip.emplace_back(2 * j + 1, 2 * j, -dRw_dpj);
					trip.emplace_back(2 * j + 1, upS, -dRw_dS);
				}
			}
			else if (wj && wj->type == Well::Injector) {
				if (!wc || wc->type != Well::Injector) {
					trip.emplace_back(2 * c + 1, 2 * c, dRw_dpi);
					trip.emplace_back(2 * c + 1, 2 * j, dRw_dpj);
					trip.emplace_back(2 * c + 1, upS, dRw_dS);
				}
			}
			else {
				if (!skipS_c) {
					trip.emplace_back(2 * c + 1, 2 * c, dRw_dpi);
					trip.emplace_back(2 * c + 1, 2 * j, dRw_dpj);
					trip.emplace_back(2 * c + 1, upS, dRw_dS);
				}
				if (!skipS_j) {
					trip.emplace_back(2 * j + 1, 2 * c, -dRw_dpi);
					trip.emplace_back(2 * j + 1, 2 * j, -dRw_dpj);
					trip.emplace_back(2 * j + 1, upS, -dRw_dS);
				}
			}
		}
	}
	Eigen::SparseMatrix<double> J(2 * N, 2 * N);
	J.setFromTriplets(trip.begin(), trip.end());
	return J;
}


inline NewtonReport newtonSolve(const Model& m, Eigen::VectorXd& x,
	const std::vector<double>& Sw_old, double dt, int max_iter = 25, double tol = 1e-6) {

	Eigen::VectorXd R = assembleResidual(m, x, Sw_old, dt);
	double r0 = R.norm();
	if (r0 == 0.0) return { true, 0 };

	for (int it = 1; it <= max_iter; ++it) {
		Eigen::SparseMatrix<double> J = assembleJacobian(m, x, dt);
		Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
		lu.compute(J);
		if (lu.info() != Eigen::Success) return { false, it };
		Eigen::VectorXd dx = lu.solve(-R);

		double f0 = 0.5 * R.squaredNorm();
		double alpha = 1.0;
		const double c1 = 1e-4;
		const double beta = 0.5;
		Eigen::VectorXd x_new = x;
		Eigen::VectorXd R_new = R;
		for (int ls = 0;ls < 20;++ls) {
			x_new = x + alpha * dx;
			double Swc = m.fluid.kr.Swc, Sor = m.fluid.kr.Sor;
			for (int c = 0; c < m.grid.numCells(); ++c) {
				const Well* w = m.wellAt(c);
				if (w && w->type == Well::Injector) continue;
				x_new(2 * c + 1) = std::min(1.0 - Sor, std::max(Swc, x_new(2 * c + 1)));
			}
			R_new = assembleResidual(m, x_new, Sw_old, dt);
			if (0.5 * R_new.squaredNorm() <= (1.0 - c1 * alpha) * f0) break;
			alpha *= beta;
		}
		x = x_new; R = R_new;
		int nCells = static_cast<int>(R.size() / 2);
		auto pIdx = Eigen::seq(0, 2 * nCells - 2, 2);
		auto sIdx = Eigen::seq(1, 2 * nCells - 1, 2);
		Eigen::VectorXd R_p = R(pIdx);
		Eigen::VectorXd R_s = R(sIdx);
		double p_res = R_p.norm();
		double s_res = R_s.norm();
		//std::cout << "P_res " << p_res << " S_res " << s_res << std::endl;
		if (R.norm() < tol * r0 + tol) {
			return { true, it };
		}
	}
	return { false, max_iter };
}