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
}}
