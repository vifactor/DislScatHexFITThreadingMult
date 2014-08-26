/*
 * ANASampleCub.cpp
 *
 *  Created on: 12 july 2014
 *      Author: Viktor Kopp
 */

#include "ANASampleCub.h"

ANASampleCub::ANASampleCub(double thickness, double size)
{
	m_thickness = thickness; m_size = size;
}

ANASampleCub::~ANASampleCub()
{
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		delete m_layers[i];
	}
}

double ANASampleCub::T_threading(double r, double phi) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		result += m_layers[i]->T(r, phi);
	}
	return result;
}

void ANASampleCub::addThreadingLayer(double rho, double b_edge, double b_screw, double rc,
		double Qx, double Qz, double nu)
{
    m_layers.push_back(new ANAThreadingLayerCub(rho, b_edge, b_screw, rc, Qx, Qz, nu));
}

void ANASampleCub::resetThreadingLayer(size_t i, double rho, double rc)
{
	if(i < m_layers.size())
	{
		m_layers[i]->m_rho = rho;
		m_layers[i]->m_Rc = rc;
	}
}
