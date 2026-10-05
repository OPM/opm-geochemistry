/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OPM_GEOCHEM_INTERFACE_HPP
#define OPM_GEOCHEM_INTERFACE_HPP

#include <opm/simulators/geochemistry/GeoChemIF.h>

#include <nlohmann/json.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


class OpmGeoChemInterface final : public CGeoChemIF
{
public:
    /*!
    * \brief Constructor
    */
    OpmGeoChemInterface() : CGeoChemIF()
    { }

    /*!
    * \brief Initialize geochemistry solver from standard input file
    *
    * \param file_name Input file name
    * \param temperature Temperature
    * \param porosity Porosity
    * \param user_order Species order
    *
    * \note user_order ensures that internal geochemistry solver and OPM orders species concentration the same way
    */
    void initialize(const std::string& file_name,
                    double temperature,
                    double porosity,
                    const std::vector<std::string>& user_order = {});

    /*!
    * \brief Initialize geochemistry solver from JSON file
    *
    * \param file_name JSON file name
    * \param user_order Species order
    *
    * \note user_order ensures that internal geochemistry solver and OPM orders species concentration the same way
    */
    void initialize_json(const std::string& file_name,
                         const std::vector<std::string>& user_order = {});

    /*!
    * \brief Initialize geochemistry solver from JSON file and species from OPM deck
    *
    * \param file_name JSON file name
    * \param species Species names
    * \param minerals Minerals names
    * \param ion_ex Ion exchange species names
    * \param charge_balance Require charge balance in equilibrium solver
    * \param tol Tolerances for material blance and pH
    * \param splay_tree_resolution Splay tree resolution
    */
    void initialize_from_opm_deck(const std::string& file_name,
                                  const std::vector<std::string>& species,
                                  const std::vector<std::string>& minerals,
                                  const std::vector<std::string>& ion_ex,
                                  bool charge_balance,
                                  std::pair<double, double> tol,
                                  int splay_tree_resolution);

    /*!
    * \brief Map minerals in the OPM deck to the order of the geochemistry solver
    *
    * \param deck_minerals Mineral names in deck order
    * \return Index in the geochemistry solver's mineral order of each mineral in deck_minerals
    *
    * \warning Must be called after initialize_from_opm_deck().
    */
    std::vector<std::size_t> mineral_index_map(const std::vector<std::string>& deck_minerals) const;

    /*!
    * \brief Calculate initial mineral concentration
    *
    * \param Cmin Output vector for mineral concentration
    * \param porosity Porosity of the rock
    * \param weight_mineral Weight percent of mineral in rock
    */
    void calculate_initial_mineral_concentration(std::vector<double>& Cmin,
                                                 double porosity,
                                                 const std::unordered_map<std::string, double> &weight_mineral);

    /*!
    * \brief Set surface concentrations
    *
    * \param swat Water saturation
    * \param C_tot Concentration vector
    * \param frac_DL Diffusion layer fraction
    * \param C_io Ion exchange concentration
    *
    * Surface concentrations calculated and added to the total concentration vector. Uses the underlying SetSurfaceConc()
    * function.
    */
    void set_surface_concentrations(double swat,
                                    std::vector<double>& C_tot,
                                    double& frac_DL,
                                    const std::unordered_map<std::string, double>& C_io);

    /*!
    * \brief Get the log a mineral object
    *
    * \return Reference to log a mineral vector
    */
    std::vector<double>& get_log_a_mineral();
    const std::vector<double>& get_log_a_mineral() const;

private:
    /*!
    * \brief Append species from OPM deck to geochemistry JSON file
    *
    * \param file_name JSON file name
    * \param charge_balance Bool for requiring charge balance on H+/pH
    * \param species Species names
    * \param minerals Mineral names
    * \param ion_ex Ion exchange names
    * \returns JSON object with definitions for geochemistry solver
    */
    static nlohmann::json opmDeckSpeciesToJSON_(const std::string& file_name,
                                                bool charge_balance,
                                                const std::vector<std::string>& species = {},
                                                const std::vector<std::string>& minerals = {},
                                                const std::vector<std::string>& ion_ex = {});

    /*!
    * \brief Check if species order are the same as internally in geochemistry solver
    *
    * \param file_name JSON file name
    * \param user_order Species order
    */
    void checkUserOrder_(const std::vector<std::string>& user_order,
                         const std::optional<std::string>& file_name = std::nullopt);

    /*!
    * \brief Set tolerances for solver
    *
    * \param GCS_in Reference to chemistry solver
    * \param mbal_tol Material balance tolerance
    * \param ph_tol pH tolerance
    */
    void set_solver_tolerances(GCSolver& GCS_in,
                               double mbal_tol,
                               double ph_tol);
};

#endif //OPM_GEOCHEM_INTERFACE_HPP
