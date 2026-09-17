#pragma once

// equation between nucleation rate and estimated temperature.

// function:  T2 = 1 / (1/T1 - ln(r2/r1) * R/Ea)
// inverse:  r2 = r1 * exp((Ea/R) * (1/T1 - 1/T2)) 
// r1, r2 are transition counts per TRANSITION_INTERVAL.
// T1, T2 are in K
// Ea: the activation enegergy
// R: gas constant.

#include <cmath>

namespace nucleation
{
	constexpr double GAS_CONSTANT = 8.31;

	// Count per interval needed at temperature t2, given count r1 at temperature t1.
	inline double rate_from_temperature(double t2, double r1, double t1, double ea)
	{
		return r1 * std::exp((ea / GAS_CONSTANT) * (1.0 / t1 - 1.0 / t2));
	}

	inline double temperature_from_rate(double r2, double r1, double t1, double ea)
	{
		if (r2 <= 0 || r1 <= 0) return 0;
		double inv_t = 1.0 / t1 - std::log(r2 / r1) * GAS_CONSTANT / ea;
		if (inv_t <= 0) return 0;
		return 1.0 / inv_t;
	}
}
