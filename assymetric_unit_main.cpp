#include "../../CASMcode/include/casm/external/Eigen/Core"
#include "../../CASMcode/include/casm/external/Eigen/Dense"
#include <algorithm>
#include <casmutils/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include "assymetric_unit.hpp"
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
	
    CLI::App app{"This utility reads in a POSCAR and prints to screen a list of the symmetrically unique atoms in the structure"};
	casmutils::fs::path structurepath;
	CLI::Option* structure_path=app.add_option("-s, --structure", structurepath, "Path to structure in POSCAR format, required")-> required();
	
	CLI11_PARSE(app, argc, argv);
	
    
    std::cout<<"The chosen POSCAR is"<<structurepath<< std::endl;
	casmutils::xtal::Structure input_structure=casmutils::xtal::Structure::from_poscar(structurepath);

    std::vector<casmutils::sym::CartOp> factor_group=make_factor_group(input_structure);
    std::vector<Eigen::Vector3d> assymmetric_unit=make_assymmetric_unit(input_structure.basis_sites(), factor_group,tol);
    std::cout<<"Assymmetric Unit Sites are located at:"<<std::endl;
    for(int i=0; i<assymmetric_unit.size();i++){
        std::cout<<"Site "<<i<<":  "<<assymmetric_unit[i]<<std::endl;
    }
	return 0; 
}

