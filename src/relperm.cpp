#pragma once
#include <algorithm>
#include <cmath>

struct BrooksCorey {

	double krw0 = 1.0;
	double kro0 = 1.0;
	double Swc = 0.2;
	double Sor = 0.2;
	double n = 2.0;

	double Se(double Sw) const {
		double s = (Sw - Swc) / (1.0 - Swc - Sor);
		return std::min(1.0, std::max(0.0, s));
	}

	double dSe_dSw() const {
		return 1.0 / (1.0 - Swc - Sor);
	}

	double krw(double Sw) const {
		double s = Se(Sw);
		return krw0 * std::pow(s, n);
	}

	double kro(double Sw) const {
		double s = Se(Sw);
		return kro0 * std::pow(1.0 - s, n);
	}

	double dkrw_dSw(double Sw) const {
		double s = Se(Sw);
		if (s <= 0.0 || s >= 1.0) return 0.0;
		return krw0 * n * std::pow(s, n - 1.0) * dSe_dSw();
	}

	double dkro_dSw(double Sw) const {
		double s = Se(Sw);
		if (s <= 0.0 || s >= 1.0) return 0.0;
		return -kro0 * n * std::pow(1.0 - s, n - 1.0) * dSe_dSw();
	}


};