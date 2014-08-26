/*
 * ANASampleCub.h
 *
 *  Created on: 12 july 2014
 *      Author: kopp
 */

#ifndef ANASAMPLECUB_H_
#define ANASAMPLECUB_H_

#include "ANAThreadingLayerCub.h"
#include <vector>

class ANASampleCub
{
public:
	ANASampleCub(double thickness, double size);
	~ANASampleCub();

	void addThreadingLayer(double rho, double b_edge, double b_screw, double rc,
			double Qx, double Qz, double nu);
	void resetThreadingLayer(size_t i, double rho, double rc);
	double T_threading(double r, double phi) const;
	
	size_t getNbThreadingLayers() const {return m_layers.size();}
	double m_thickness;
	double m_size;
protected:
	std::vector<ANAThreadingLayerCub * > m_layers;
};

#endif /* ANASampleCub_H_ */
