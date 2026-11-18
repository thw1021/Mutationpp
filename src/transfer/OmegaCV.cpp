/**
 * @file OmegaCV.cpp
 *
 * @brief Implementation of OmegaCV.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include "Mixture.h"
#include "TransferModel.h"

namespace Mutation {
    namespace Transfer {

/**
* Represents a coupling between chemistry and vibrational energy modes.
*/
class OmegaCV : public TransferModel
{
/**
* The necessary functions for the computation of the  source terms for
* Vibration-Chemistry-Vibration coupling are implemented.
*
* The dissociation of molecules is a procedure which removes energy from the vibrational mode
* of them and re-introduces it as chemical energy of the produced atoms and molecules. In order to accurately determine
* the exact energy lost from the vibrational mode in a dissociation reaction one must use state-to-state models in order to
* account for the fact that not all the vibrational states have the exact same probability to dissociate. Generally, the higher
* vibrational levels have a greater probability to dissociate compared to the lower ones because the latter ones should ladder climb
* before dissociating. Nonetheless, the state-to-state models are not easy to use in a CFD code due to the increase of the computational
* cost. Therefore, one should rely on accurate Boltzmann equilibrium models in order to take into account for this mode of energy transfer.
*
* Taking into account the above, the models for the Vibration-Chemistry-Vibration energy transfer can be generaly placed in to two different categories,
* the non-preferential and the preferential models. In the first category, all the energy states have the same probability to dissociate compared to the
* preferential ones which assume higher probability of dissociation for the higher energy levels. Currently, only the simplest non-preferential model has been
* implemented, but the code gives the opportunity to be extended for more complicated Vibration-Chemistry-Vibration energy transfer models.
*
*/

public:
	OmegaCV(Mutation::Mixture& mix)
		: TransferModel(mix)
	{
		m_ns = mixture().nSpecies();
		mp_wrk1 = new double [m_ns];
		mp_wrk2 = new double [m_ns];
		mp_wrk3 = new double [m_ns];
		mp_dT   = new double [m_ns];
		mp_dTv  = new double [m_ns];
		mp_cpv  = new double [m_ns];
		mp_cpel = new double [m_ns];
		mp_jac  = new double [m_ns*m_ns];
	};

	~OmegaCV()
	{
		delete [] mp_wrk1;
		delete [] mp_wrk2;
		delete [] mp_wrk3;
		delete [] mp_dT;
		delete [] mp_dTv;
		delete [] mp_cpv;
		delete [] mp_cpel;
		delete [] mp_jac;
	};
/**
 * Computes the source terms of the Vibration-Chemistry energy transfer in \f$ [J/(m^3\cdot s)] \f$
 *
 * \f[ \Omega^{CV}_{mi} = \sum c_1 e^V_{mi} \dot{\omega}_i\f]
 *
 */
	double source()
	{
		static int i_transfer_model = 0;
		switch (i_transfer_model){
		   case 0:
			  return compute_source_Candler();
		  break;
		   default:
			  std::cerr << "The selected Chemistry-Vibration-Chemistry model is not implemented yet";
			  return 0.0;
		}
	}

private:
	int m_ns;
	double* mp_wrk1;
	double* mp_wrk2;
	double* mp_wrk3;
	double* mp_dT;
	double* mp_dTv;
	double* mp_cpv;
	double* mp_cpel;
	double* mp_jac;

	double const compute_source_Candler();

public:
	/**
	 * Analytic derivatives of Omega_CV =
	 * Ru*T*sum_i hv_i(Tv) wdot_i(T,rho)/M_i with respect to T, Tv and rho_s.
	 *
	 * d/dT  = Ru*sum hv_i wdot_i/M_i + Ru*T*sum hv_i (dwdot_i/dT)/M_i
	 * d/dTv = Ru*sum cpv_i wdot_i/M_i     (hv_i = Hv_i/(RuT), dHv_i/dTv = cpv_i*Ru)
	 * d/dr = Ru*T*sum hv_i (dwdot_i/drho_s)/M_i
	 * The chemistry derivatives use the analytic rate-law path (no finite
	 * differences), so this returns the exact derivative of this model.
	 */
	virtual bool sourceJacobian(double& dT, double& dTv, double* const dRho)
	{
		const double T = mixture().T();

		mixture().speciesHOverRT(NULL, NULL, NULL, mp_wrk1, NULL, NULL);
		mixture().netProductionRates(mp_wrk2);
		mixture().netProductionRatesJacobian(mp_dT, mp_dTv);
		mixture().jacobianRho(mp_jac);
		mixture().speciesCpOverR(
			T, mixture().Te(), mixture().Tr(), mixture().Tv(),
			mixture().Tel(), NULL, NULL, NULL, mp_cpv, mp_cpel);

		dT = 0.0;
		dTv = 0.0;
		for (int s = 0; s < m_ns; ++s)
			dRho[s] = 0.0;

		for (int i = 0; i < m_ns; ++i) {
			const double mw = mixture().speciesMw(i);
			const double hv = mp_wrk1[i];
			const double wd = mp_wrk2[i];

			// hv_i = Hv_i/(Ru*T): the explicit T of the prefactor and the 1/T
			// inside hv_i cancel, so Omega_CV = sum Hv_i wdot_i/M_i and only
			// the chemistry derivative survives in d/dT.
			dT  += RU*T*(hv*mp_dT[i]/mw);
			dTv += RU*(mp_cpv[i]*wd/mw);

			for (int s = 0; s < m_ns; ++s)
				dRho[s] += RU*T*(hv*mp_jac[i*m_ns+s]/mw);
		}
		return true;
	}
};

 /**
 * Non-preferential Model according to Candler with
 *
 * with \f$ c_1 \f$ equal to 1 for non-preferential models. For more information:
 *
 * G. V. Candler, MacCormack, Computation of weakly ionized hypersonic flows in thermochemical nonequilibrium, Journal of
 * Thermophysics and Heat Transfer, 1991, 5(11):266
 *
 */

double const OmegaCV::compute_source_Candler()
{


	 // Getting Vibrational Energy
	 mixture().speciesHOverRT(NULL, NULL, NULL, mp_wrk1, NULL, NULL);

	 // Getting Production Rate
	 mixture().netProductionRates(mp_wrk2);

	 // Inner Product
	 double c1 = 1.0E0;
	 double sum = 0.E0;

	 for(int i = 0 ; i < m_ns; ++i)
		 sum += mp_wrk1[i]*mp_wrk2[i]/mixture().speciesMw(i);

	 return(c1*sum*mixture().T()*RU);
 }

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaCV, TransferModel> omegaCV("OmegaCV");
      
    } // namespace Transfer
} // namespace Mutation
