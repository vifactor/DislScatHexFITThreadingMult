/*
 * ANAThreadingLayerHex.h
 *
 *  Created on: 26 aug. 2014
 *      Author: kopp
 */

#ifndef ANATHREADINGLAYERHEX_H_
#define ANATHREADINGLAYERHEX_H_

#include <gsl/gsl_math.h>
#include <iostream>

class ANAThreadingLayerHex
{
public:
	ANAThreadingLayerHex(double rho, double b_edge, double b_screw, double rc, double Qx, double Qz, double nu);
	virtual ~ANAThreadingLayerHex();

	double T(double r, double phi) const;
	double m_rho;
	double m_Rc;
protected:
	double m_gb_screw, m_gb2_screw;
	double m_gb_edge, m_gb2_edge;

	double m_C_screw, m_C1_edge, m_C2_edge;

	double T_screw(double r) const;
	double T_edge(double r, double phi) const;
	inline double chi_screw() const;
	inline double chi_edge(double phi) const;
};

#endif /* ANAThreadingLayerHex_H_ */
