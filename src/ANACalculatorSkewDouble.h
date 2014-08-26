/*
 * ANACalculatorSkewDouble.h
 *
 *  Created on: 24 july 2014
 *      Author: kopp
 */

#ifndef ANACALCULATORSKEWDOUBLE_H_
#define ANACALCULATORSKEWDOUBLE_H_

#include "ANASampleCub.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>

class ANACalculatorSkewDouble
{
public:
	ANACalculatorSkewDouble(double Qx, double Qz, 
	        double lambda = 0.154, // Cu K-alpha wavelength in nm
	        double epsabs = 1e-5 // absolute error for integration
	        );
	virtual ~ANACalculatorSkewDouble();

	virtual void setResolution(double fwhm_qx, double fwhm_qz);
	virtual void setSample(ANASampleCub * );
	ANASampleCub * getSample() {return m_sample;}
	virtual double I(const double omega) const;
	virtual double I(const double omega, double & error) const;
	
	friend double ana_skew_double_integrand(double x, void *params);
protected:
    ANASampleCub * m_sample;

    double m_alpha, m_angcoef;
	double m_resol2_x, m_resol2_z;
	mutable double m_frequency;

	inline double T_threading(double x) const;


	/*gsl integration staff*/
	gsl_integration_workspace * m_workspace;
	gsl_integration_workspace * m_cyclic_workspace;
	gsl_integration_qawo_table * m_qawo_table;
	mutable gsl_function m_function;
	size_t m_limit;
	double m_epsabs;
};

#endif /* ANACALCULATORSKEWDOUBLE_H_ */
