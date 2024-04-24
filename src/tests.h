#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quadrature(int order, bool border)
        {
            double w = 0;
            Quadrature q;
            q = q.get_quadrature(order, border);
            std::cout<<"test quadrature"<<"\n";
            std::cout<<"nb_points" << q.nb_points() << "\n";
            for (int pt = 0; pt < q.nb_points(); pt++){
            	w += q.weight(pt);
            	std::cout<< "x : " << q.point(pt).x << " et y : " << q.point(pt).y << "\n";
            	std::cout<< "w : " << w <<"\n";
            	}
            std::cout<<w<<"\n";
            return true;
        
        }
        
        bool test_map(std::string M, bool border, int i)
        {
	// lire square.mesh et créer le maillage comme dans les test précédents 
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping triangle(mesh, border, i);
            return true;
        }
        
        bool test_transform()
        {
             Mesh mesh;
             mesh.load("data/square.mesh");
             ElementMapping triangle(mesh, false, 4);
             
             vertex point;
             point.x = 0.2;
             point.y = 0.4;
             
	     vertex R = triangle.transform(point);
	     std::cout << "après transormation, x = " << R.x << " et y = "<< R.y << "\n";
             return true;
        }
        
        bool test_Jacobian_Matrix()
        {
             Mesh mesh;
             mesh.load("data/square.mesh");
             ElementMapping triangle(mesh, false, 4);
             
             vertex point;
             point.x = 0.2;
             point.y = 0.4;
             
             DenseMatrix Jout = triangle.jacobian_matrix(point);
             Jout.print();
             return true;
    	}
    	
    	bool test_Jacobian_Det()
    	{
    	     Mesh mesh;
             mesh.load("data/square.mesh");
             ElementMapping triangle(mesh, false, 4);
             
             vertex point;
             point.x = 0.2;
             point.y = 0.4;
             
             double det;
             det = triangle.jacobian(point);
             std::cout << "le déterminant de la matrice est : " << det << ". \n";
             return true;
        }
        
        bool test_ShapeFunction (int dim, int order)
        {
             ShapeFunctions SF(dim, order);
             return true;
        }
        
        bool test_nb_functions()
        {
             ShapeFunctions SF1(2,1);
             std::cout << "il s'agit d'un triangle " << SF1.nb_functions()  << "\n";
             ShapeFunctions SF2(1,1);
             std::cout << "il s'agit d'un segment " << SF2.nb_functions() << "\n";
             return true;
        }
        
        bool test_evaluate (int i)
        {
       	     ShapeFunctions SF(2,1);
       	     vertex point; 
       	     point.x = 0.2;
             point.y = 0.4;
             
             std::cout << "la valeur de phi au point : " << i << " vaut : " << SF.evaluate(i, point) << "\n";
             
             return true;
        }
        
        bool test_G_evaluate(int i)
        {
             ShapeFunctions SF(2,1);
             vertex point;
             point.x = 0.2;
             point.y = 0.4;
        	
             std::cout << "le gradient de phi au point " << i << " est x = " << SF.evaluate_grad(i, point).x << " y = " << SF.evaluate_grad(i, point).y << " \n";
             return true;
        }
        
        double unit_fct( vertex v )
        {
            return 1.;
        }
        
        bool test_AEM()
        {
        	Mesh mesh;
             	mesh.load("data/square_fine.mesh");
            	ElementMapping elt_mapping(mesh, false, 4);
        	
        	ShapeFunctions reference_functions (2,1);
        	
        	Quadrature quadrature;
        	quadrature  = quadrature.get_quadrature(0, false);
        	vertex pt_quad = quadrature.point(0);

        	DenseMatrix Ke;
		Ke.set_size(reference_functions.nb_functions(),reference_functions.nb_functions());
	
        	
        	assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, unit_fct, Ke);
        	Ke.print();
		return true;
	}
        
        bool test_LtGMatrix()
        {// Necessaire pour obtenir Ke
        	Mesh mesh; // maillage
        	mesh.load("data/square_fine.mesh");
        	ElementMapping elt_mapping(mesh, false, 4);
        	ShapeFunctions reference_functions (2,1);
        	Quadrature quadrature;
        	quadrature  = quadrature.get_quadrature(0, false);
        	vertex pt_quad = quadrature.point(0);
        	DenseMatrix Ke;
		Ke.set_size(reference_functions.nb_functions(),reference_functions.nb_functions());
		assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, unit_fct, Ke);
        	Ke.print();
        	
	//Construction de la matrice K
		
		int tailleK = mesh.nb_vertices();
		SparseMatrix K(tailleK);
		int t = 4; // index du triangle
		local_to_global_matrix(mesh, t, Ke, K);
		K.print();
		return true;
	}
		
	bool test_aDBc()
	{//Préparation pour appliquer les conditions de D
		Mesh mesh; // maillage
        	mesh.load("data/square_fine.mesh");
        	ElementMapping elt_mapping(mesh, false, 4);
        	ShapeFunctions reference_functions (2,1);
        	Quadrature quadrature;
        	quadrature  = quadrature.get_quadrature(0, false);
        	vertex pt_quad = quadrature.point(0);
        	DenseMatrix Ke;
		Ke.set_size(reference_functions.nb_functions(),reference_functions.nb_functions());
		assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, unit_fct, Ke);
		int tailleK = mesh.nb_vertices();
		SparseMatrix K(tailleK);
		int t = 4; // index du triangle
		local_to_global_matrix(mesh, t, Ke, K);

		//initialisation attribute_is_dirichlet
		std::vector<bool> attribute_is_dirichlet;
		attribute_is_dirichlet.push_back(true);
		
		//Attribution d'un attribut à tout les edges
		mesh.set_attribute(unit_fct, 0, true);
		
		//initialisation values
		std::vector<double> values;
		//values[mesh.nb_vertices()];
		for (int i = 0; i < mesh.nb_vertices(); i++)
		{
			//values[i] = mesh.get_vertex(i).x + mesh.get_vertex(i).y;
			values.push_back(mesh.get_vertex(i).x + mesh.get_vertex(i).y);
			//std::cout<<values[i]<< "\n";
		}
		
		
		//initialiation F
		std::vector<double> F(mesh.nb_vertices(), 0);
		
		apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F);
		for (int k = 0; k<mesh.nb_vertices(); k++)
		{
			std::cout<<F[k] << "\n";
		}
		return true;
	}
        	
}}
