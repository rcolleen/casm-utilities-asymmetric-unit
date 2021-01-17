#include "assymetric_unit.hpp"
#include <algorithm>
#include <cmath>


std::vector<Eigen::Vector3d> make_orbit(const Eigen::Vector3d& coordinate,
					const std::vector<casmutils::sym::CartOp>& symops,  
					const casmutils::xtal::Lattice& lattice)
{
	double tol=1e-5;
	std::vector<Eigen::Vector3d> orbit;
        for (const casmutils::sym::CartOp& symop: symops)
	{
            Eigen::Vector3d transformedcoord = symop * coordinate;
            VectorPeriodicCompare_f make_sure_no_repeated_coords(transformedcoord, tol, lattice);
	    if (std::find_if(orbit.begin(), orbit.end(), make_sure_no_repeated_coords)==orbit.end())
	//    //Need condition to check if the created coord is not exactly the same as another coord VectorPeriodicCompare_f maybe
			    {
			    orbit.emplace_back(transformedcoord);	
			    }    
        }
	return orbit;
}

// Indexes vectors of the different asymmetric orbits of the system. This cam further be broken into the asymmetric atoms in the system
// For reference these use the interstitial coordinates after the old ones are taken out
std::vector<int>
label_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const std::vector<casmutils::sym::CartOp>& symops,
				  const casmutils::xtal::Lattice& lattice,
                                  double tol)
{
    int label=0;
    std::vector<int> coordinate_tags(coordinates.size(),-1);
    for(int ix = 0; ix<(coordinate_tags.size()); ++ix)
    {
        if(coordinate_tags[ix]!=-1)
        {
            continue;
        }

        const Eigen::Vector3d& coord=coordinates[ix];
        std::vector<Eigen::Vector3d> coord_orbit=make_orbit(coord, symops,lattice);

        for(int ixx=ix; ixx<coordinate_tags.size(); ++ixx)
        {
            if(coordinate_tags[ixx]!=-1)
            {
                continue;
            }

            VectorPeriodicCompare_f equals_ixx_coord(coordinates[ixx],tol,lattice);
	    if(std::find_if(coord_orbit.begin(),coord_orbit.end(),equals_ixx_coord)!=coord_orbit.end())
            {
                coordinate_tags[ixx]=label;	
            }
        }
	++label;
    }

    for(auto l : coordinate_tags)
    {
        assert(l!=-1);
    }
	
    return coordinate_tags;
}

std::vector<std::vector<Eigen::Vector3d>>
bin_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const std::vector<casmutils::sym::CartOp>& symops,
		                  const casmutils::xtal::Lattice& lattice,
                                  double tol)
{
	std::vector<int> coordinate_tags= label_by_symmetrical_equivalence(coordinates, symops, lattice, tol);
	//goes through length of coordinate list and connects coordinate to orbit	
	int num_orbits=*std::max_element(coordinate_tags.begin(), coordinate_tags.end())+1;
	std::vector<std::vector<Eigen::Vector3d>> orbit_container;
	orbit_container.resize(num_orbits);
       	for (int i=0; i<coordinates.size(); i++)
	{
		Eigen::Vector3d temp_coord=coordinates[i];
		int label=coordinate_tags[i];
		orbit_container[label].emplace_back(temp_coord);
	}
	return orbit_container;	
}

std::vector<Eigen::Vector3d> make_asymmetric_unit(const std::vector<Eigen::Vector3d>& complete_structure_basis,
                                                  const std::vector<casmutils::sym::CartOp>& symops,
						  const casmutils::xtal::Lattice& lattice,
                                                  double tol)
{
	std::vector<Eigen::Vector3d> asymmetric_unit_collated;
	std::vector<std::vector<Eigen::Vector3d>> total_orbits=bin_by_symmetrical_equivalence(complete_structure_basis, symops, lattice, tol);
	for (const auto& orbit: total_orbits)
	{
	if (orbit.size()>0)
	  {  
			asymmetric_unit_collated.push_back(orbit[0]);
	  }
	else
		std::cout<<"orbit of size zero here!"<<std::endl;	
	}
	return asymmetric_unit_collated;
}

