/**
 * @file RateManager.h
 *
 * @brief Declaration of RateManager class.
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

#ifndef KINETICS_RATEMANAGER_H
#define KINETICS_RATEMANAGER_H

#include "RateLawGroup.h"
#include "Reaction.h"

namespace Mutation {
    namespace Kinetics {

class Reaction;

/**
 * Manages the efficient computation of rate coefficients for the evaluation of 
 * reaction rates.
 */
class RateManager
{
public:

    /**
     * Creates a new RateManager which manages the reaction rate coefficients
     * for the given reactions.
     */
    RateManager(size_t ns, const std::vector<Reaction>& reactions);
    
    /**
     * Destructor.
     */
    ~RateManager();
    
    /**
     * Updates the current values of the rate coefficients.
     */
    void update(const Thermodynamics::Thermodynamics& thermo);

    /**
     * Updates the temperature derivatives of the rate coefficients (logarithmic
     * derivatives with respect to T and Tv) at the current state.  Intended to
     * be called immediately before the derivative accessors below; less than
     * repeated finite differences by construction.
     */
    void updateDerivative(const Thermodynamics::Thermodynamics& thermo);

    /// \f$ \partial \ln k_f / \partial T \f$ at fixed species densities.
    const double* const dlnkf_dT() const { return mp_dlnkf_dT; }
    /// \f$ \partial \ln k_f / \partial T_v \f$ at fixed species densities.
    const double* const dlnkf_dTv() const { return mp_dlnkf_dTv; }
    /// \f$ \partial \ln k_b / \partial T \f$ at fixed species densities.
    const double* const dlnkb_dT() const { return mp_dlnkb_dT; }
    /// \f$ \partial \ln k_b / \partial T_v \f$ at fixed species densities.
    const double* const dlnkb_dTv() const { return mp_dlnkb_dTv; }
    
    /**
     * Returns a pointer to the forward rate coefficients evaluated at the 
     * forward temperature.
     */
    const double* const lnkf() { return mp_lnkf; }
    
    /**
     * Returns a pointer to the forward rate coefficients evaluated at the 
     * forward temperature.
     */
    const double* const lnkb() { return mp_lnkb; }
    
    /**
     * Returns the indices of irreversible reactions.
     */
    const std::vector<size_t>& irrReactions() const {
        return m_irr;
    }

private:

    /**
     * Add a new reaction whose rate law should be managed.  This method will
     * select the appropriate temperature to evaluate the given reaction's rate
     * law based on the type of reaction it is.
     */
    void addReaction(const size_t rxn, const Reaction& reactions);

    /**
     * Helper function which enumerates all of the possible reaction types and
     * calls the appropriate addRate() template based on the reaction type.
     */
    template <int NReactionTypes>
    void selectRate(const size_t rxn, const Reaction& reaction);

    /**
     * Small helper function which adds a new rate law to the manager based on
     * which rate law groups the forward and reverse rates fall into.
     */
    template <typename ForwardGroup, typename ReverseGroup>
    void addRate(const size_t rxn, const Reaction& reaction);

private:

    /// Number of species in the mixture
    const size_t m_ns;
    
    /// Number of reactions in the mechanism
    const size_t m_nr;

    /// Collection of different rate law groups (can be forward or reverse)
    RateLawGroupCollection m_rate_groups;
    
    /// Storage for forward rate coefficient at forward temperature
    double* mp_lnkf;
    
    /// Storage for forward rate coefficient at backward temperature
    double* mp_lnkb;
    
    /// Storage for species Gibbs free energies
    double* mp_gibbs;

    /// Temperature-derivative block (11*nr doubles):
    /// [0,nr) dlnkf/dT, [nr,2nr) dlnkf/dTv, [2nr,3nr) dlnkb/dT,
    /// [3nr,4nr) dlnkb/dTv, [4nr,5nr) dlnkeq/dTg,
    /// [5nr,7nr) dlnkf/dTg (forward and reverse groups),
    /// [7nr,9nr) chain dTg/dT, [9nr,11nr) chain dTg/dTv
    double* mp_deriv_block;
    double* mp_dlnkf_dT;
    double* mp_dlnkf_dTv;
    double* mp_dlnkb_dT;
    double* mp_dlnkb_dTv;
    double* mp_dlnkeq_dTg;
    double* mp_dlnkf_dTg;
    double* mp_chain_dT;
    double* mp_chain_dTv;
    
    /// Stores the indices for which the forward and reverse temperature
    /// evaluations are equal
    std::vector<size_t> m_to_copy;
    
    /// Stores the indices of non-reversible reactions
    std::vector<size_t> m_irr;
};


    } // namespace Kinetics
} // namespace Mutation


#endif // KINETICS_RATEMANAGER_H
