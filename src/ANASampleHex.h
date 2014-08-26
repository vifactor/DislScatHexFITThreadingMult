/*
 * ANASampleHex.h
 *
 *  Created on: 12 july 2014
 *      Author: kopp
 */

#ifndef ANASAMPLEHEX_H_
#define ANASAMPLEHEX_H_

#include "ANAThreadingLayerHex.h"
#include <vector>

class ANASampleHex
{
public:
	ANASampleHex(double thickness, double size);
	~ANASampleHex();

	void addThreadingLayer(double rho, double b_edge, double b_screw, double rc,
			double Qx, double Qz, double nu);
	void resetThreadingLayer(size_t i, double rho, double rc);
	double T_threading(double r, double phi) const;
	
	size_t getNbThreadingLayers() const {return m_layers.size();}
	double m_thickness;
	double m_size;
protected:
	std::vector<ANAThreadingLayerHex * > m_layers;
};

#endif /* ANASampleHex_H_ */
