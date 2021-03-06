#include <algorithm>
#include <casmutils/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include "asymmetric_unit.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <CLI/CLI.hpp>
#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"


//user inputted mesh function declaration

int main(int argc, char* argv[]) { 
	double tol=1E-6;
    bool cart_flag;
	
    CLI::App app{"This utility reads in a POSCAR and prints to screen a list of the symmetrically unique atoms in the structure"};
	casmutils::fs::path structurepath;
	CLI::Option* structure_path=app.add_option("-s, --structure", structurepath, "Path to structure in POSCAR format, required")-> required();
    app.add_flag("-c, --cart", cart_flag, "Cartesian coordinates out, default is fractional");	
	CLI11_PARSE(app, argc, argv);
	
    
    std::cout<<"The chosen POSCAR is "<<structurepath<< std::endl;
	casmutils::xtal::Structure input_structure=casmutils::xtal::Structure::from_poscar(structurepath);

    std::vector<casmutils::sym::CartOp> factor_group=make_factor_group(input_structure, tol);
    std::vector<Eigen::Vector3d> basis;
    for (int i=0; i<input_structure.basis_sites().size(); i++)
    {
        //generate list of basis in cartersian coordinates for asym 
        basis.push_back(input_structure.basis_sites()[i].cart());
    }
    std::vector<Eigen::Vector3d> asymmetric_unit=make_asymmetric_unit(basis, factor_group, input_structure.lattice(), tol);
    std::cout<<"Assymmetric Unit Sites are located at:  "<<std::endl;

    if(cart_flag){
        for(int i=0; i<asymmetric_unit.size();i++){
        std::cout<<"Site "<<i<<":  "<<asymmetric_unit[i].transpose()<<std::endl;
        }
    }
    else{
        for(int i=0; i<asymmetric_unit.size();i++){
            Eigen::Vector3d coordinate = casmutils::xtal::cartesian_to_fractional(asymmetric_unit[i], input_structure.lattice());
            
        std::cout<<"Site "<<i<<":  "<<coordinate.transpose()<<std::endl;
        }
    }
	return 0; 
}

