#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    
    const bool t_quadrature=false;
    const bool t_map = false;
    const bool t_transform = false;
    const bool t_Mjacobienne = false;
    const bool t_DetJacobienne = false;
    
    const bool t_SF = false;
    const bool t_nb_functions = false;
    const bool t_evaluate = false;
    const bool t_G_evaluate = false;
    const bool t_AEM = false;
    const bool t_LtGMatrix = true;

    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    
    if (t_quadrature)  Tests::test_quadrature(0, false);
    if (t_map) Tests::test_map("data/square.mesh", false, 4);
    if (t_transform) Tests::test_transform();
    if (t_Mjacobienne) Tests::test_Jacobian_Matrix();
    if (t_DetJacobienne) Tests::test_Jacobian_Det();
    
    if (t_SF) Tests::test_ShapeFunction(1,1);
    //if (t_SF) Tests::test_ShapeFunction(3,1);
    if (t_nb_functions) Tests::test_nb_functions();
    if (t_evaluate) Tests::test_evaluate(2);
    if (t_G_evaluate) Tests::test_G_evaluate(0);
    
    if (t_AEM) Tests::test_AEM();
    if (t_LtGMatrix) Tests::test_LtGMatrix();
}

void run_simu()
{

    const bool simu_pure_dirichlet = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) {
        Simu::pure_dirichlet_pb("data/square.mesh", verbose);
    }
}

int main( int argc, const char * argv[] )
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
