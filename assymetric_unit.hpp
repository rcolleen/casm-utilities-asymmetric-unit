#ifndef INTERSTITIAL_MESH_H
#define INTERSTITIAL_MESH_H

#include "casmutils/xtal/lattice.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/symmetry.hpp"
#include "casmutils/xtal/structure.hpp"
#include <vector>
#include "vectorfunctions.hpp"

//This makes an orbit (or a set of symmetrically equivalent atoms)
std::vector<Eigen::Vector3d> make_orbit(const Eigen::Vector3d& coordinates, const std::vector<casmutils::sym::CartOp>& symops, const casmutils::xtal::Lattice& lattice);


//This function labels atoms by their respective orbits
std::vector<int>
label_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const std::vector<casmutils::sym::CartOp>& symops,
                                  const casmutils::xtal::Lattice& lattice,
                                  double tol);

//This function is used to bin sets of atoms into different orbits
std::vector<std::vector<Eigen::Vector3d>>
bin_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const std::vector<casmutils::sym::CartOp>& symops,
                                  const casmutils::xtal::Lattice& lattice,
                                  double tol);

//finds the asymmetric unit of a ser of atoms through application of the factor group
std::vector<Eigen::Vector3d> make_asymmetric_unit(const std::vector<Eigen::Vector3d>& complete_structure_basis, 
                                  const std::vector<casmutils::sym::CartOp>& symops,
				  const casmutils::xtal::Lattice& lattice, 
			          double tol);
#endif 
