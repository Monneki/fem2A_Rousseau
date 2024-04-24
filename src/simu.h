#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        	{
            	//std::cout << "Solving a pure Dirichlet problem" << std::endl;
            	if ( verbose ) 
            		{
                	std::cout << " Simulation de la solution d'un pb de Dirichlet pure sur le maillage square_fine. \n" << "Pour modifier le maillaige en question il faut modifier dans main.cpp le fichier prit en entier, et dans simu.h le nom du fichier en sortie.\n " << "!! il faut que le nom du fichier en sortie .bb soit le même que celui prit en entrée " << "\n";
            		}
            	Mesh mesh;
            	mesh.load(mesh_filename);
        	ShapeFunctions reference_functions(2,1);
        	Quadrature quadrature;
        	quadrature = quadrature.get_quadrature(0,false);
        	//initialisation matrice K
        	int tailleK = mesh.nb_vertices();
		SparseMatrix K(tailleK);
		std::vector<double> F(mesh.nb_vertices(), 0);
		
		
		
		
		for (int triangle = 0; triangle<mesh.nb_triangles(); triangle ++)
		{
			ElementMapping elt_mapping(mesh, false, triangle);
			//Initialisation Ke
			DenseMatrix Ke;
        		Ke.set_size(reference_functions.nb_functions(),reference_functions.nb_functions());
			assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, unit_fct, Ke);
			local_to_global_matrix(mesh, triangle, Ke, K);
		}
		//initialisation condition Dirichlet
		std::vector<bool> attribute_is_dirichlet;
		attribute_is_dirichlet.push_back(true);
		mesh.set_attribute(unit_fct, 0, true);
		//initialisation values
		std::vector<double> values;
		for (int i = 0; i < mesh.nb_vertices(); i++)
		{
			values.push_back(mesh.get_vertex(i).x + mesh.get_vertex(i).y);
		}
		apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
		std::vector<double> X(mesh.nb_vertices(), 0);
		std:solve(K, F, X);
		std::string export_name = "data/square_fine";
		mesh.save(export_name+".mesh");
		save_solution(X, export_name+".bb");

        }

    }

}
