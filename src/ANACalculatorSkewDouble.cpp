/*
 * ANACalculatorSkewDouble.cpp
 *
 *  Created on: 24 july 2014
 *      Author: kopp
 */

#include "ANACalculatorSkewDouble.h"

double ana_skew_double_integrand(double x, void *params)
{
	ANACalculatorSkewDouble * calculator;
	static double result;

	calculator = static_cast<ANACalculatorSkewDouble *> (params);

	result = exp(-calculator->T_threading(x))
			* cos(calculator->m_frequency * x);

	return result;
}

ANACalculatorSkewDouble::ANACalculatorSkewDouble(double Qx, double Qz, 
                                                double lambda, double epsabs)
{
    double Q;
    double sinThetaB, cosThetaB,
            sinPsi, cosPsi, 
            sinPhi, cosPhi,
            cosAlpha;
    
	m_sample = NULL;

	m_resol2_x = 0.0;
	m_resol2_z = 0.0;

	m_frequency = 0.0;
	
	Q = sqrt(Qx * Qx + Qz * Qz);
	/*skew geometry angles*/
	sinThetaB = Q * lambda / (4 * M_PI);
	cosThetaB = sqrt(1.0 - gsl_pow_2(sinThetaB));
	sinPsi = Qz / Q;
	sinPhi = sinPsi * sinThetaB;
	cosPsi = sqrt(1.0 - gsl_pow_2(sinPsi));
	cosPhi = sqrt(1.0 - gsl_pow_2(sinPhi));
	cosAlpha = sinThetaB * cosPsi / cosPhi;
	
	m_angcoef = Q * cosThetaB / cosPhi;
    m_alpha = acos(cosAlpha);

	/*initialization of integration staff*/
	m_limit = 10000;
	m_epsabs = epsabs;
	m_function.function = &ana_skew_double_integrand;
	m_function.params = this;

	m_qawo_table = gsl_integration_qawo_table_alloc(0.0, 0.0, GSL_INTEG_COSINE,
			m_limit);
	m_workspace = gsl_integration_workspace_alloc(m_limit);
	m_cyclic_workspace = gsl_integration_workspace_alloc(m_limit);

    /*switch off exception like behavior on bad integral convergence*/
	gsl_set_error_handler_off ();
}

ANACalculatorSkewDouble::~ANACalculatorSkewDouble()
{
	if(m_qawo_table)
		gsl_integration_qawo_table_free(m_qawo_table);
	if(m_workspace)
		gsl_integration_workspace_free(m_workspace);
	if(m_cyclic_workspace)
		gsl_integration_workspace_free(m_cyclic_workspace);
}

void ANACalculatorSkewDouble::setSample(ANASampleCub * sample)
{
    m_sample = sample;
}

void ANACalculatorSkewDouble::setResolution(double fwhm_qx, double fwhm_qz)
{
	/* FWHM = 2 * sqrt(2 * log(2)) / sigma */
	/* resol2 = 1 / (2 * sigma**2) */
	m_resol2_x = gsl_pow_2(fwhm_qx / 4) / log(2.0);
	m_resol2_z = gsl_pow_2(fwhm_qz / 4) / log(2.0);
}

double ANACalculatorSkewDouble::T_threading(double x) const
{
	static double result;

	result = m_sample->T_threading(fabs(x), m_alpha);

	return result;
}

double
ANACalculatorSkewDouble::I(const double omega, double & error) const
{
	static double result, abserr;

	m_frequency = m_angcoef * omega;

	//gsl_integration_qawo_table_set(m_qawo_table, m_frequency, 1, GSL_INTEG_COSINE);
	//gsl_integration_qawf (&m_function, m_a, m_epsabs, m_limit, m_workspace, m_cyclic_workspace, m_qawo_table, &result, &abserr);

	gsl_integration_qagiu (&m_function, 0.0, m_epsabs, 1.e-6, m_limit, m_workspace, &result, &abserr);
	result *= 2;
	error = 2 * abserr;

	return result;
}

double ANACalculatorSkewDouble::I(const double omega) const
{
	static double result, abserr;

	result = I(omega, abserr);

	return result;
}
