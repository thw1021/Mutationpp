/**
 * @file RateLawGroup.h
 *
 * Defines the abstract class RateLawGroup which is used by the RateManager
 * class to combine rate laws which are evaluated at the same temperature in
 * order to compute the rate coefficient for these reactions efficiently.
 *
 * @see class RateLawGroup
 * @see class RateLawGroup1T
 * @see class RateLawGroupCollection
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

#ifndef KINETICS_RATE_LAW_GROUP_H
#define KINETICS_RATE_LAW_GROUP_H

#include <map>
#include <typeinfo>
#include <vector>

#include "RateLaws.h"
#include "Reaction.h"
//#include "StateModel.h"
#include "StoichiometryManager.h"

class StateModel;

namespace Mutation {
    namespace Kinetics {

/**
 * Abstract base class which defines the interface for all RateLawGroup objects
 * which evaluate like rate raws to increase efficiency.
 */
class RateLawGroup
{
public:

    /**
     * Constructor.
     */
    RateLawGroup() : m_last_t(-1.0) {}

    /**
     * Destructor.
     */
    virtual ~RateLawGroup() { };

    /**
     * Adds a new rate to evaluate with this group.
     */
    virtual void addRateCoefficient(
        const size_t rxn, const RateLaw* const p_rate) = 0;
    
    /**
     * Adds a reaction which uses the backward temperature represented by this
     * group.
     */
    void addReaction(const size_t rxn, const Reaction& reaction) {
        m_reacs.addReaction(rxn, reaction.reactants());
        m_prods.addReaction(rxn, reaction.products());
    }
    
    /**
     * Returns the temperature used in the last evaluation of the rate
     * coefficients.
     */
    double getT() const { return m_t; }
    
    /**
     * Evaluates all of the rates in the group and stores in the given vector.
     */
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk) = 0;
        
    /**
     * Computes \Delta G / RT for this rate law group and subtracts these values
     * for each of the reactions in this group.
     */
    void subtractLnKeq(size_t ns, double* const p_g, double* const p_r) const
    {        
        // Compute G_i/RT - ln(Patm/RT)
        const double val = std::log(ONEATM / (RU * m_t));
        for (int i = 0; i < ns; ++i)
            p_g[i] -= val;
        
        // Now subtract \Delta[G_i/RT - ln(Patm/RT)]_j
        m_reacs.decrReactions(p_g, p_r);
        m_prods.incrReactions(p_g, p_r);
    }

    /**
     * Accumulates the group-temperature derivative of the equilibrium-constant
     * term that converts k_f into k_b, for every reaction whose reverse rate
     * coefficient is evaluated in this group.
     *
     * ln k_b = ln k_f + \Delta G/RT - \Delta\nu \ln(P_{atm}/R_uT) with
     * \partial (G_i/R_uT)/\partial T = -H_i/(R_uT^2) and
     * \partial \ln T/\partial T = 1/T, so the contribution is
     * \f$ \sum_{react} (H_i/R_uT - 1)/T - \sum_{prod} (H_i/R_uT - 1)/T \f$
     * evaluated at the group temperature.
     *
     * @note requires m_t to be set by a preceding lnk()/dlnkdT() call.
     */
    void accumulateDlnKeq(
        const size_t ns, const Thermodynamics::Thermodynamics& thermo,
        double* const p_dlnkeq) const
    {
        // D75 (SU2 NEMO): differentiate the same G/RT that subtractLnKeq()
        // evaluates, instead of using the analytic -H/(R T^2).  The electronic
        // Boltzmann factors are tabulated (RrhoDB m_use_tables) and the table
        // interpolates Z and T^2 dZ/dT independently, so -H/(R T^2) is not the
        // derivative of the tabulated G/RT (relative error up to 2e-4 at
        // 8000-12000 K).  See SU2 Docs/NEMOReport/reports/
        // NMP-D75_keq_table_rootcause_2026-11-23.md.
        std::vector<double> mu(ns), gp(ns), gm(ns);
        const double hstep = 1.0e-5 * m_t;
        thermo.speciesSTGOverRT(m_t + hstep, gp.data());
        thermo.speciesSTGOverRT(m_t - hstep, gm.data());
        for (size_t i = 0; i < ns; ++i)
            mu[i] = -(gp[i] - gm[i]) / (2.0 * hstep) + 1.0 / m_t;

        m_reacs.incrReactions(mu.data(), p_dlnkeq);
        m_prods.decrReactions(mu.data(), p_dlnkeq);
    }

    /**
     * Evaluates the group-temperature derivative of the rate coefficients and
     * of the k_f -> k_b conversion, and returns the chain factors that convert
     * a group-temperature derivative into a derivative with respect to T and Tv.
     *
     * Arrays are indexed by the same reaction index used in lnk() (a reversed
     * rate coefficient of reaction j is stored at j+nReactions), except
     * p_dlnkeq_dTg which is indexed by the reaction index j and is additive.
     */
    virtual void dlnkdT(
        const Thermodynamics::Thermodynamics& thermo,
        const Thermodynamics::StateModel* const p_state,
        double* const p_dlnkf_dTg, double* const p_chain_dT,
        double* const p_chain_dTv, double* const p_dlnkeq_dTg) = 0;


protected:

    /// This is the temperature computed to evaluate the rate law (should be set
    /// in the lnk() function)
    double m_t;
    double m_last_t;
    
    /// Stores the reactants for reactions that will use this rate law for the
    /// reverse direction
    StoichiometryManager m_reacs;
    
    /// Stores the products for reactions that will use this rate law for the
    /// reverse direction
    StoichiometryManager m_prods;
};


/**
 * Groups reaction rate laws based on a single temperature that are the same
 * kind and evaluated at the same temperature together so that they may be 
 * evaluated efficiently.
 */
template <typename RateLawType, typename TSelectorType>
class RateLawGroup1T : public RateLawGroup
{
public:

    /**
     * Adds a new rate to evaluate with this group.
     */
    virtual void addRateCoefficient(
        const size_t rxn, const RateLaw* const p_rate)
    {
        m_rates.push_back(
            std::make_pair(rxn, dynamic_cast<const RateLawType&>(*p_rate))
        );
    }

    /**
     * Evaluates all of the rates in the group and stores in the given vector.
     */
    virtual void lnk(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk)
    {
        // Determine the reaction temperature for this group
        m_t = TSelectorType().getT(p_state);

        // Update only if the temperature has changed
        //if (std::abs(m_t - m_last_t) > 1.0e-10) {
            const double lnT  = std::log(m_t);
            const double invT = 1.0 / m_t;

            for (int i = 0; i < m_rates.size(); ++i) {
                const std::pair<size_t, RateLawType>& rate = m_rates[i];
                p_lnk[rate.first] = rate.second.getLnRate(lnT, invT);
            }
        //}

        // Save this temperature
        m_last_t = m_t;
    }

    /**
     * Evaluates dln(k)/dTg for every rate in this group and writes the chain
     * factors dTg/dT and dTg/dTv, plus the additive equilibrium-constant term
     * for the reverse rates handled by this group.
     */
    virtual void dlnkdT(
        const Thermodynamics::Thermodynamics& thermo,
        const Thermodynamics::StateModel* const p_state,
        double* const p_dlnkf_dTg, double* const p_chain_dT,
        double* const p_chain_dTv, double* const p_dlnkeq_dTg)
    {
        m_t = TSelectorType().getT(p_state);

        const double lnT  = std::log(m_t);
        const double invT = 1.0 / m_t;
        const double dTdT  = TSelectorType().dTdT(p_state);
        const double dTdTv = TSelectorType().dTdTv(p_state);

        /* D72: in a single-temperature model Tv and Te are not independent
         * variables - the state model pins them to T (Tv() == Te() == T() by
         * construction), so the rate temperature T_g = T_g(T, Tv, Te) responds
         * to a change in T through *every* selector channel, while it has no
         * independent Tv channel at all (a finite difference that perturbs Tv
         * alone leaves the state model unchanged, hence dT_g/dTv = 0).
         * Without this collapse ParkSelector returned dT_g/dT = 0.5 instead of
         * 1 and dT_g/dTv = 0.5 instead of 0, and TeSelector returned
         * dT_g/dT = 0 instead of 1, so every Park and electron-temperature
         * reaction carried a wrong temperature derivative (and a spurious
         * d(omega)/dTv). */
        double chainT  = dTdT;
        double chainTv = dTdTv;
        if (thermo.nEnergyEqns() == 1) {
            chainT  = dTdT + dTdTv;
            chainTv = 0.0;
        }

        for (int i = 0; i < m_rates.size(); ++i) {
            const std::pair<size_t, RateLawType>& rate = m_rates[i];
            p_dlnkf_dTg[rate.first] = rate.second.dlnkdT(lnT, invT);
            p_chain_dT[rate.first]  = chainT;
            p_chain_dTv[rate.first] = chainTv;
        }

        accumulateDlnKeq(thermo.nSpecies(), thermo, p_dlnkeq_dTg);

        m_last_t = m_t;
    }

private:

    /// vector of rates to evaluate
    std::vector< std::pair<size_t, RateLawType> > m_rates;
};


/**
 * Small helper class which provides comparison between two std::type_info
 * pointers.
 */
struct CompareTypeInfo {
    bool operator ()(const std::type_info* a, const std::type_info* b) const {
        return a->before(*b);
    }
};


/**
 * Manages a collection of RateLawGroup objects such that only one object of any
 * RateLawGroup type is ever created in each collection.
 */
class RateLawGroupCollection
{
public:

    typedef std::map<const std::type_info*, RateLawGroup*, CompareTypeInfo>
        GroupMap;

    /**
     * Destructor.
     */
    ~RateLawGroupCollection()
    {
        GroupMap::iterator iter = m_group_map.begin();
        for ( ; iter != m_group_map.end(); ++iter) {
            delete iter->second;
            iter->second = NULL;
        }
    }
    
    /**
     * Returns the number of different rate law groups in this collection.
     */
    size_t nGroups() const { return m_group_map.size(); }
    
    /**
     * Returns the GroupMap managed by this RateLawGroupCollection object.
     */
    const GroupMap& groups() const { return m_group_map; }

    /**
     * Adds a new rate law to be managed by this collection of rate law groups.
     */
    template <typename GroupType>
    void addRateCoefficient(const size_t rxn, const RateLaw* const p_rate)
    {
        if (m_group_map[&typeid(GroupType)] == NULL)
            m_group_map[&typeid(GroupType)] = new GroupType();
        m_group_map[&typeid(GroupType)]->addRateCoefficient(rxn, p_rate);
    }
    
    /**
     * Adds a reaction to the manager which allows for the calculation of the 
     * \Delta G / RT term for reverse rate coefficients.
     */
    template <typename GroupType>
    void addReaction(const size_t rxn, const Reaction& reaction)
    {
        if (m_group_map[&typeid(GroupType)] == NULL)
            m_group_map[&typeid(GroupType)] = new GroupType();
        m_group_map[&typeid(GroupType)]->addReaction(rxn, reaction);
    }

    /**
     * Computes the rate coefficients in this collection and stores them in
     * the vector at the index corresponding to their respective reaction.
     */
    void logOfRateCoefficients(
        const Thermodynamics::StateModel* const p_state, double* const p_lnk)
    {
        // Compute the forward rate constants
        GroupMap::iterator iter = m_group_map.begin();
        for ( ; iter != m_group_map.end(); ++iter)
            iter->second->lnk(p_state, p_lnk);
    }

    /**
     * Evaluates dln(k)/dTg and the chain factors for every group in the
     * collection.  p_dlnkeq_dTg must be zeroed by the caller and is indexed by
     * the reaction index.
     */
    void dlnkdT(
        const Thermodynamics::Thermodynamics& thermo,
        const Thermodynamics::StateModel* const p_state,
        double* const p_dlnkf_dTg, double* const p_chain_dT,
        double* const p_chain_dTv, double* const p_dlnkeq_dTg)
    {
        GroupMap::iterator iter = m_group_map.begin();
        for ( ; iter != m_group_map.end(); ++iter)
            iter->second->dlnkdT(
                thermo, p_state, p_dlnkf_dTg, p_chain_dT, p_chain_dTv,
                p_dlnkeq_dTg);
    }
    
    /**
     * Subtracts ln(keq) from the provided rate coefficients.
     */
    void subtractLnKeq(
        const Thermodynamics::Thermodynamics& thermo, double* const p_g,
        double* const p_lnk)
    {
        const size_t ns = thermo.nSpecies();
        GroupMap::iterator iter = m_group_map.begin();
        for ( ; iter != m_group_map.end(); ++iter) {
            const RateLawGroup* p_group = iter->second;
            thermo.speciesSTGOverRT(p_group->getT(), p_g);
            p_group->subtractLnKeq(ns, p_g, p_lnk);
        }
    }

private:
    
    /// Collection of RateLawGroup objects
    GroupMap m_group_map;
};

    } // namespace Kinetics
} // namespace Mutation

#endif // KINETICS_RATE_LAW_GROUP_H

