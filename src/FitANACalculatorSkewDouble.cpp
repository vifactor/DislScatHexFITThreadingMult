/*
 * FitANACalculatorSkewDouble.cpp
 *
 *  Created on: 25 july 2013
 *      Author: kopp
 */

#include "FitANACalculatorSkewDouble.h"
using namespace boost;

FitANACalculatorSkewDouble::FitANACalculatorSkewDouble(
                    const std::vector<ANACalculatorSkewDouble *>& calculators)
{
	m_calculators = calculators;
	m_scales.resize(m_calculators.size());
    m_backgrounds.resize(m_calculators.size());
    
    /*all samples contain the same number of layers and interfaces*/
    size_t nblayers = m_calculators[0]->getSample()->getNbThreadingLayers();
    initParameterNames(nblayers);
}

void FitANACalculatorSkewDouble::initParameterNames(size_t nblayers)
{
    /*Data fit parameter names*/
    m_scale_names.resize(m_calculators.size());
    m_background_names.resize(m_calculators.size());
    for(size_t id = 0; id < m_calculators.size(); ++id)
    {
        m_scale_names.at(id) = "Data.[" + lexical_cast<std::string>(id) + "].I0";
        m_background_names.at(id) = "Data.[" + lexical_cast<std::string>(id) + "].Ibg";
    }
    /*theading layers names*/
    m_th_density_names.resize(nblayers);
    m_th_rc_names.resize(nblayers);
    for(size_t id = 0; id < nblayers; ++id)
    {
        m_th_density_names.at(id) = "Sample.dislocations.threading.[" 
                            + lexical_cast<std::string>(id) 
                            + "].rho";
        m_th_rc_names.at(id) = "Sample.dislocations.threading.[" 
                            + lexical_cast<std::string>(id) 
                            + "].rc";
    }
}

void FitANACalculatorSkewDouble::reinit(const NonlinearFit::CalculatorParameterMap& params)
{
    static double rho_mf;
	static double rho_th, rc_th;
    
    /*update each calculator*/
    for(size_t cid = 0; cid < m_calculators.size(); ++cid)
    {
        /*reinitialization scale and background coefficients*/
        m_scales[cid] = params.find(m_scale_names[cid])->second;
        m_backgrounds[cid] = params.find(m_background_names[cid])->second;


        /*
         * reinitialization of densities and correlation radii of threading dislocations.
         * initially threading dislocation density is given in [cm-2]
         * coefficient 1e-14 transforms it to [nm-2]
         */
        for(size_t lid = 0;
                lid < m_calculators[cid]->getSample()->getNbThreadingLayers();
                ++lid)
        {
            rho_th = params.find(m_th_density_names[lid])->second * 1e-14;
            rc_th = params.find(m_th_rc_names[lid])->second;
            m_calculators[cid]->getSample()->resetThreadingLayer(lid, rho_th, rc_th);
            
            
            //DEBUG
            /*
            std::cout << "calculator:\t" << cid << std::endl;
            std::cout << m_th_density_names[lid] << ":\t" << rho_th << std::endl;
            std::cout << m_th_rc_names[lid]  << ":\t" << rc_th << std::endl;
            */
        }
	}
}

double FitANACalculatorSkewDouble::eval(const NonlinearFit::CalculatorArgument * arg)
{
	static double omega, result;
	static int id;

	omega = static_cast<const ANACalculatorSkewDoubleArgument* >(arg)->m_omega;
	id = static_cast<const ANACalculatorSkewDoubleArgument* >(arg)->m_id;
    
	result = m_scales[id] * m_calculators[id]->I(omega) + m_backgrounds[id];

	return result;
}

