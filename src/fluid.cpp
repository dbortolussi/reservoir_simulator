#pragma once
#include "relperm.cpp"

struct Fluid {
	double mu_w = 1.0e-3;
	double mu_o = 5.0e-3;
	BrooksCorey kr;

	double lambda_w(double Sw) const {
		return kr.krw(Sw) / mu_w;
	}

	double lambda_o(double Sw) const {
		return kr.kro(Sw) / mu_o;
	}

	double lambda_t(double Sw) const {
		return lambda_w(Sw) + lambda_o(Sw);
	}

	double fw(double Sw) const {
		double lw = lambda_w(Sw);
		double lt = lambda_t(Sw);
		return (lt > 0.0) ? lw / lt : 0.0;
	}

	double dlambda_t_dSw(double Sw) const {
		return kr.dkrw_dSw(Sw) / mu_w + kr.dkro_dSw(Sw) / mu_o;
	}

	double dfw_dSw(double Sw) const {
		double lw = lambda_w(Sw);
		double lt = lambda_t(Sw);

		double dlw = kr.dkrw_dSw(Sw) / mu_w;
		double dlt = dlambda_t_dSw(Sw);

		return (lt > 0.0) ? (dlw*lt-lw*dlt)/(lt*lt) : 0.0;
	}
};