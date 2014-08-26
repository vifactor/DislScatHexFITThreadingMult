/*
 * FitANACalculatorSkewDouble.h
 *
 *  Created on: 24 july 2014
 *      Author: kopp
 */

#ifndef FITANACALCULATORSKEWDOUBLE_H_
#define FITANACALCULATORSKEWDOUBLE_H_

#include "ANACalculatorSkewDouble.h"
#include "NonlinearFit.h"
#include <boost/lexical_cast.hpp>

class ANACalculatorSkewDoubleArgument : public NonlinearFit::CalculatorArgument
{
public:
	ANACalculatorSkewDoubleArgument(double omega, int id)
	{
	    m_omega = omega;
	    m_id = id;
	}
	virtual ~ANACalculatorSkewDoubleArgument() {}
	
	double m_omega;
	int m_id;
};

class FitANACalculatorSkewDouble:
		public NonlinearFit::FitCalculator
{
public:
	FitANACalculatorSkewDouble(
	        const std::vector<ANACalculatorSkewDouble *>& calculators);
	virtual ~FitANACalculatorSkewDouble() {}
	virtual void reinit(const NonlinearFit::CalculatorParameterMap& params);
	virtual double eval(const NonlinearFit::CalculatorArgument * arg);
protected:
	std::vector<ANACalculatorSkewDouble *> m_calculators;
	std::vector<double> m_scales, m_backgrounds;
	std::vector<std::string> m_scale_names,
	                         m_background_names,
	                         m_th_density_names,
	                         m_th_rc_names;
	
	void initParameterNames(size_t nblayers);
};

#endif /* FITANACALCULATORSKEWDOUBLE_H_ */
