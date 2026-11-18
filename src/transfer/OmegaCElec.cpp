/*
 * @file OmegaCElec.cpp
 *
 * @brief Implementation of OmegaCElec.
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
 * Represents a coupling between chemistry and electronic energy modes.
 *
 * \f[ \sum_{i \in \varepsilon } \dot{\omega}_i e^E_i \f].
 *
 * where the index \f$ \varepsilon \f$ represents the set oh heavy species
 *
 * More information about the above model can be found in @cite Panesi 2008
 */
class OmegaCElec : public TransferModel
{
public:
	OmegaCElec(Mutation::Mixture& mix)
		: TransferModel(mix)
	{
		m_ns = mix.nSpecies();
		mp_wrk1 = new double [m_ns];
		mp_wrk2 = new double [m_ns];
		mp_dT   = new double [m_ns];
		mp_dTv  = new double [m_ns];
		mp_cpv  = new double [m_ns];
		mp_cpel = new double [m_ns];
		mp_jac  = new double [m_ns*m_ns];
	};

	~OmegaCElec()
	{
		delete [] mp_wrk1;
		delete [] mp_wrk2;
		delete [] mp_dT;
		delete [] mp_dTv;
		delete [] mp_cpv;
		delete [] mp_cpel;
		delete [] mp_jac;
	};

	double source()
	{
		mixture().speciesHOverRT(NULL, NULL, NULL, NULL, mp_wrk1, NULL);
		mixture().netProductionRates(mp_wrk2);

		double sum = 0.0;

		for(int i = 0 ; i < mixture().nSpecies(); ++i)
			sum += mp_wrk1[i]*mp_wrk2[i]/mixture().speciesMw(i);

		return (sum*mixture().T()*RU);
	}

private:
	int m_ns;
	double* mp_wrk1;
	double* mp_wrk2;
	double* mp_dT;
	double* mp_dTv;
	double* mp_cpv;
	double* mp_cpel;
	double* mp_jac;

public:
	/**
	 * Analytic derivatives of Omega_CElec =
	 * Ru*T*sum_i hel_i(Tel) wdot_i(T,rho)/M_i.
	 * Same structure as OmegaCV with the electronic enthalpy/cp.
	 */
	virtual bool sourceJacobian(double& dT, double& dTv, double* const dRho)
	{
		const double T = mixture().T();

		mixture().speciesHOverRT(NULL, NULL, NULL, NULL, mp_wrk1, NULL);
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
			const double hel = mp_wrk1[i];
			const double wd  = mp_wrk2[i];

			// hel_i = Hel_i/(Ru*T): the prefactor T and the 1/T inside hel_i
			// cancel, leaving only the chemistry derivative in d/dT.
			dT  += RU*T*(hel*mp_dT[i]/mw);
			dTv += RU*(mp_cpel[i]*wd/mw);

			for (int s = 0; s < m_ns; ++s)
				dRho[s] += RU*T*(hel*mp_jac[i*m_ns+s]/mw);
		}
		return true;
	}
};

// Register the transfer model
Mutation::Utilities::Config::ObjectProvider<
    OmegaCElec, TransferModel> omegaCElec("OmegaCElec");

    } // namespace Transfer
} // namespace Mutation
