/*
 *
 *  Created on: 30 may 2012
 *  Modified on: 25 july 2014
 *      Author: Viktor Kopp
 */

#include "Engine.h"

using namespace NonlinearFit;
using namespace Geometry;
using namespace boost;

Vector3d toThreadingFrame(const Vector3d& burgers,
                     const Vector3d& normal)
{
    Vector3d norm_dir;
    
    norm_dir = normalize(normal);
    
    /*screw and edge components of the dislocation*/
    double b_screw = inner_prod(burgers, norm_dir);
    double b_edge = norm_2(burgers - b_screw * norm_dir);
    
    return Vector3d(b_edge, 0.0, b_screw);
}

Vector3d toCoplanarFrame(const Vector3d& Q, const Vector3d& normal)
{
    Vector3d norm_dir;
    double Qx, Qy, Qz;

    norm_dir = normalize(normal);
    
    Qz = inner_prod(Q, norm_dir);
    Qy = 0.0;//in coplanar geometry Qy = 0
    Qx = norm_2(Q - Qz * norm_dir);
    
    return Vector3d(Qx, Qy, Qz);
}

/*----------- end handy functions -----------*/
Engine::Engine()
{
	m_fitter = NULL;
	m_programSettings = NULL;
	m_fit_calculator = NULL;
}

Engine::~Engine()
{
	if(m_programSettings)
		delete m_programSettings;

    init();
}

void Engine::exec(const filesystem::path & workDir)
{
    m_WorkDir = workDir;
	m_programSettings = new ProgramSettings;
	m_programSettings->read(m_WorkDir);

	init();
	setupComponents();
	doWork();
	saveResult();
}

void Engine::init()
{
	for(size_t i = 0; i < m_calculators.size(); ++i)
	{
        if(m_calculators[i])
        {
            if(m_calculators[i]->getSample())
                delete m_calculators[i]->getSample();
            delete m_calculators[i];
        }
	}
	if(m_fit_calculator)
		delete m_fit_calculator;
	if(m_fitter)
		delete m_fitter;

	/*vector of fit parameters with their final values*/
	m_fParametersFinal.clear();

	m_exp_intens_vals.clear(); m_ini_intens_vals.clear(); m_fin_intens_vals.clear();
	m_domega_vals.clear();

	for(size_t i = 0; i < m_DataPoints.size(); ++i)
	{
		delete m_DataPoints[i].m_Argument;
	}
	m_DataPoints.clear();
	
	m_vals_sizes.clear();
	m_vals_sizes.push_back(0);
}

void Engine::doWork()
{
    calcI(m_ini_intens_vals);
    /*perform fit and file updates only if nbinterations > 0 */
    if(m_programSettings->getFitConfig().nbIterations > 0)
    {
        fitI();
        //update configuration files
        updateConfigFile(m_programSettings->getSampleConfigFile());
        updateConfigFile(m_programSettings->getDataConfigFile());
        //save resume file
        saveResume();
        //calculate final intensities
        calcI(m_fin_intens_vals);
    }
    else
        m_fin_intens_vals = m_ini_intens_vals;
}

void Engine::calcI(std::vector<double>& vals)
{
	m_fit_calculator->reinit(m_programSettings->getCPMap());

	for(size_t ipt = 0; ipt < m_DataPoints.size(); ipt++)
		vals[ipt] = m_fit_calculator->eval(m_DataPoints[ipt].m_Argument);
}

filesystem::path
Engine::getOutFilename(const filesystem::path& infile) const
{
    filesystem::path outfile;
    
    outfile = infile.filename();
	outfile.replace_extension("ft");
	return m_WorkDir / outfile;
}

void Engine::saveResult() const
{
    for(size_t id = 0; id < m_programSettings->getDataConfig().size(); ++id)
        saveFitData(id);
}

void Engine::saveFitData(size_t id) const
{
    filesystem::path filename;
    
    filename = getOutFilename(m_programSettings->getDataConfig(id).file);
    
    std::ofstream fout(filename.c_str());
	fout << "#qx\tqz\tIexp\tIini\tIfin" << std::endl;
	for(size_t ipt = m_vals_sizes[id]; ipt < m_vals_sizes[id + 1]; ipt++)
	{
        fout << m_domega_vals[ipt]  << "\t"
			<< m_exp_intens_vals[ipt] << "\t"
			<< m_ini_intens_vals[ipt] << "\t"
			<< m_fin_intens_vals[ipt] << std::endl;
	}
	fout.close();
}

void Engine::fitI()
{
    //carry out fit
    //with max "calculationSettings.modefitSettings.nbIterations" iterations
    m_fitter->fit(NonlinearFitter::fitLIN,
        m_programSettings->getFitConfig().nbIterations);
    printFitterInfo(m_fitter);

    //get new values (best fit values) of fit parameters
    m_fParametersFinal = m_fitter->getFitParameters();

    /*update cpMap parameters with newly fitted values*/
    mergeParameters(m_programSettings->getCPMap(),
        m_fParametersFinal);
}

void Engine::saveResume() const
{
	boost::filesystem::path filename;
	std::ofstream fout;

	/*transform rc to M = rc/rd = rc*rho^(1/2) */
	filename = m_WorkDir / m_programSettings->getResumefile();

    fout.open(filename.c_str());
    fout << "In " << m_fitter->getNbIterations() 
        << " iterations ||f||^2 reduced from "
        << m_fitter->getFinit() << " to " 
        << m_fitter->getFfin() << std::endl;
    fout << std::endl;
    for (NonlinearFit::FitParameterList::const_iterator it=m_fParametersFinal.begin();
            it!=m_fParametersFinal.end(); ++it)
    {
        fout << it->m_Name << " => " << it->m_Value << std::endl;
    }
    fout.close();
}

void Engine::updateConfigFile(const boost::filesystem::path & cfgfile) const
{
	libconfig::Config cfg;
	boost::filesystem::path oldcfgfile;

	oldcfgfile = cfgfile;
	oldcfgfile.replace_extension(".~cfg");
	/*save old configuration in a file with extension ~cfg*/
	boost::filesystem::rename(cfgfile, oldcfgfile);

	try
	{
		// Read the configuration file. If there is an error, report it
		cfg.readFile(oldcfgfile.c_str());

		//reset configuration putting the new fitted values from m_fParametersFinal
		for(size_t i = 0; i < m_fParametersFinal.size(); ++i)
		{
		    if(cfg.exists(m_fParametersFinal[i].m_Name))    
                cfg.lookup(m_fParametersFinal[i].m_Name) = 
                                                m_fParametersFinal[i].m_Value;
		}
		//save new configuration
		cfg.writeFile (cfgfile.c_str());
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Engine::Exception("I/O error while reading file:\t" + oldcfgfile.native());
	} catch (const libconfig::ParseException &pex)
	{
		throw Engine::Exception("Parse error at "
									+ lexical_cast<std::string>(pex.getFile()) 
									+ ":"
									+ lexical_cast<std::string>(pex.getLine()) 
									+ " - "
									+ lexical_cast<std::string>(pex.getError()));
	}catch(const libconfig::SettingNotFoundException &nfex)
	{
		throw Engine::Exception(nfex.getPath());
	}catch(libconfig::SettingTypeException& tex)
	{
		throw Engine::Exception(lexical_cast<std::string>(tex.getPath())
		                             + "(" 
		                             + lexical_cast<std::string>(tex.what())
		                             + ")");
	}
}

void Engine::setupComponents()
{
    for(size_t id = 0; id < m_programSettings->getDataConfig().size(); ++id)
        setupCalculator(id);
    m_fit_calculator = new FitANACalculatorSkewDouble(m_calculators);
    
    std::cout << "--- Data files ---" << std::endl;
    for(size_t id = 0; id < m_programSettings->getDataConfig().size(); ++id)
	    readData(id);
    m_ini_intens_vals.resize(m_exp_intens_vals.size(), 0.0);
    m_fin_intens_vals.resize(m_exp_intens_vals.size(), 0.0);

	setupFitter();
}

void Engine::readData(size_t id)
{
    DataReader dr("\t ");
    filesystem::path datapath;

    datapath = m_WorkDir / m_programSettings->getDataConfig(id).file;
    dr.parse(datapath.native());
    
    if(!dr.columnExists("[intensity]"))
        throw Engine::Exception("Column \"[intensity]\" has not been found in " +
                            m_programSettings->getDataConfig(id).file.native());
    else
    {
        const DataReader::ColumnType& intensity = dr.columnGet("[intensity]");
        
        if(dr.columnExists("[domega]"))
        {
            const DataReader::ColumnType& domega= dr.columnGet("[domega]");
            
            /*allocate arguments and residuals*/
            for(size_t i = 0; i < intensity.size(); ++i)
            {
                m_exp_intens_vals.push_back(intensity[i]);
                /*transform degrees to radians*/
                m_domega_vals.push_back(domega[i]);
            
                m_DataPoints.push_back(
                    NonlinearFit::DataPoint(
                            new ANACalculatorSkewDoubleArgument(domega[i] * M_PI / 180.0, id),
                            intensity[i]));
            }
        }
        else if(dr.columnExists("[omega]"))
        {
            const DataReader::ColumnType& omega= dr.columnGet("[omega]");            
            DataReader::ColumnType::const_iterator it_of_max;
            size_t index_of_max;
            double omega_max;

            //find maximum intensity == peak center
            it_of_max = std::max_element(intensity.begin(), intensity.end());
            index_of_max = it_of_max - intensity.begin();

            omega_max = omega.at(index_of_max);
            
            /*allocate arguments and residuals*/
            for(size_t i = 0; i < intensity.size(); ++i)
            {
                m_exp_intens_vals.push_back(intensity[i]);
                /*transform degrees to radians*/
                m_domega_vals.push_back(omega[i]);
            
                m_DataPoints.push_back(
                    NonlinearFit::DataPoint(
                            new ANACalculatorSkewDoubleArgument(
                                    (omega[i] - omega_max) * M_PI / 180.0,
                                    id),
                                intensity[i])
                                );
            }
        }
        else //(!dr.columnExists("[domega]") and !dr.columnExists("[omega]"))
        {
            throw Engine::Exception("Neither \"[domega]\" nor \"[omega]\" has not been found in "+
                                m_programSettings->getDataConfig(id).file.native());
        }
        
        m_vals_sizes.push_back(m_vals_sizes.back() + intensity.size());
        std::cout << "\t" << datapath << "::" << intensity.size() << std::endl;
    }
}

void Engine::setupCalculator(size_t id)
{
    Vector3d Q_vec, Q;
	Vector3d b_vec, l_vec, n_vec, b;
	MillerCubIndicesTransformator transformator(
	                                m_programSettings->getSampleConfig().a0);
	double phi;
	
	n_vec = Vector3d(0, 0, 1);

	/*transform hexagonal miller indices to vector*/
	Q_vec = transformator.toVector3d(m_programSettings->getDataConfig(id).Q);
	/*get Q in coplanar frame*/
	Q = toCoplanarFrame(Q_vec, n_vec);
	/* sign of the first index of Q_vec determines the difference between reflections
	 * like [224] and [-2-24]
	*/
	Q[0] *= GSL_SIGN (Q_vec[0]);
    
	try
	{
	    m_calculators.push_back(new ANACalculatorSkewDouble(
	            Q[0], Q[2]));
        m_calculators.back()->setSample(new ANASampleCub(
                m_programSettings->getSampleConfig().thickness,
                m_programSettings->getSampleConfig().width));
        
        m_calculators.back()->setResolution(
				m_programSettings->getDataConfig(id).resolX,
				m_programSettings->getDataConfig(id).resolZ);

		/*threading layers*/
        const std::vector<ProgramSettings::SampleConfig::ThreadingDislocationType>& 
            th_dislocations = m_programSettings->getSampleConfig().threading;
        for(size_t i = 0; i < th_dislocations.size(); ++i)
        {
            b_vec = transformator.toVector3d(th_dislocations[i].b);
            b = toThreadingFrame(b_vec, n_vec);
            m_calculators.back()->getSample()->addThreadingLayer(
                th_dislocations[i].rho * 1e-14,
                b(0), b(2),
                th_dislocations[i].rc, 
                Q(0), Q(2),
                m_programSettings->getSampleConfig().nu);
        }
	} catch (std::exception& ex)
	{
		throw Engine::Exception(
				"Calculator allocation pb (" + 
				lexical_cast<std::string>(ex.what()) + ")");
	}
}

void Engine::setupFitter()
{
	m_fitter = new PortFitter;

	m_fitter->init(m_fit_calculator,
                    m_programSettings->getFitConfig().fitParameters,
                    m_programSettings->getCPMap(),
                    m_DataPoints);
}

void Engine::printFitterInfo(NonlinearFit::NonlinearFitter * fitter)
{
	//size_t index;

	std::cout << "In " << fitter->getNbIterations() << " iterations ||f||^2 reduced from "
				<< fitter->getFinit() << " to " << fitter->getFfin() << std::endl;
	std::cout << "Nb Function evaluations:\t" << fitter->getNbFuncEval() << std::endl;
	std::cout << "Nb Jacobian evaluations:\t" << fitter->getNbGradEval() << std::endl;
	std::cout << "Reason to stop iterations:\t<" << fitter->getReasonToStop()
			<< "> Code:\t" << fitter->getReasonToStopId() << std::endl;

	std::cout << "Fitted parameters:" << std::endl;
	if(fitter->getName().compare("LEVMAR") == 0)
	{
		//index = 0;

		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "[" << i << "]\t" << fitter->getFitParameter(i).m_Name
					<< "\t" << fitter->getFitParameter(i).m_Value << " +/- "
					<< sqrt(fitter->getCovarMatrix(i, i)) << std::endl;
		}

		std::cout << "Correlation matrix:" << std::endl;
		std::cout << "\t";
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "["<< i <<"]" << "\t";
		}
		std::cout << std::endl;
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "["<< i <<"]" << "\t";
			for(size_t j = 0; j <= i; ++j)
			{
				std::cout
						<< fitter->getCovarMatrix(i, j)
								/ sqrt(
										fitter->getCovarMatrix(i, i)
												* fitter->getCovarMatrix(j, j))
						<< "\t";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		for(size_t i = 0; i < fitter->getFitParameters().size(); ++i)
		{
			std::cout << "[" << i << "]\t" << fitter->getFitParameter(i).m_Name
					<< "\t" << fitter->getFitParameter(i).m_Value << std::endl;
		}
		std::cout << "Covariance matrix has not been calculated." << std::endl;
	}
}
