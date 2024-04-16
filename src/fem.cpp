#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        //vertice_ = vecteur de vertex
        std::cout << "[ElementMapping] constructor for element " << i << " " << "\n";
        std::vector<vertex> vertices;
        if ( border ) 
        {
        	std::cout << "edge"; 
        	for (int indiceGlob = 0; indiceGlob < 2; indiceGlob++)
        	{
        		vertices_.push_back(M.get_edge_vertex(i, indiceGlob));
        	}
        	for (int v = 0; v<2; ++v)
        	{
        		std::cout<< vertices_[v].x << "\n "<< vertices_[v].y << "\n";
        	}
        }
        
        if (not border)
        {
        	std::cout << "triangle";
        	for (int indiceGlob = 0; indiceGlob < 3; indiceGlob++)
        	{
        		vertices_.push_back(M.get_triangle_vertex(i,indiceGlob));
        	}
        	/*for (int v = 0; v<3; ++v)
        	{
        		std::cout<< vertices_[v].x << "\n"<< vertices_[v].y << "\n";
        	}*/
    	}
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        std::cout << "[ElementMapping] transform reference to world space" << '\n';
        vertex r ;
        
        if (border_)
        {
        	r.x = (1-x_r.x)*vertices_[0].x + x_r.x * vertices_[1].x;
        	r.y = (1-x_r.x)*vertices_[0].y + x_r.x * vertices_[1].y;
        }
        if (not border_)
        {
        	r.x = (1 - x_r.x - x_r.y)*vertices_[0].x + x_r.x * vertices_[1].x + x_r.y * vertices_[2].x;
        	r.y = (1 - x_r.x - x_r.y)*vertices_[0].y + x_r.x * vertices_[1].y + x_r.y * vertices_[2].y;
        }
        	
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        DenseMatrix J ;
        if (border_)
        {
        	J.set_size(2,1);
        	J.set(0,0,-vertices_[0].x + vertices_[1].x);
        	J.set(1,0,-vertices_[0].y + vertices_[1].y);
        }
        if (not border_)
        {
        	J.set_size(2,2);
        	J.set(0,0,-vertices_[0].x + vertices_[1].x);
        	J.set(1,0,-vertices_[0].y + vertices_[1].y);
        	J.set(0,1,-vertices_[0].x + vertices_[2].x);
        	J.set(1,1,-vertices_[0].y + vertices_[2].y);
        }
        
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
	DenseMatrix J;
	J = jacobian_matrix(x_r);
	double det = 0;
	if (border_)
	{
		det = sqrt(J.transpose().get(0,0) * J.get(0,0) + J.transpose().get(0,1) * J.get(1,0));
		
	}
        if (not border_)
        {
        	det = J.det_2x2();
        }

        return det;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        assert(dim == 1 || dim == 2);
        assert(order == 1);
        std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
        std::cout << "[ShapeFunctions] for the order "<< order << '\n';
    }

    int ShapeFunctions::nb_functions() const
    {
        std::cout << "[ShapeFunctions] number of functions" << '\n';
        if  (dim_ == 1)
       	{
        	return 2;
       	}
        if (dim_ == 2)
        {
        	return 3;
        }
        return 0;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
	if (dim_ == 1)
	{
		double L[2] = {1 - x_r.x, x_r.x};
		return L[i];
	}
	
	if (dim_ == 2)
	{
		double L[3] =  {1 - x_r.x - x_r.y, x_r.x, x_r.y};
		return L[i];
	}
        return 0. ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
        vec2 g ;
        
        if (dim_ == 1)
        {
        	int Lxi[2] = {-1, 1};
        	int Leta[2] = { 0,0 };
        	g.x = Lxi[i];
        	g.y = Leta[i];
        }
        
        if (dim_ == 2)
        {	
        	int Lxi[3] = {-1, 1, 0};
        	int Leta[3] = {-1, 0, 1};
        	g.x = Lxi[i];
        	g.y = Leta[i];
        }
        
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )//taille Ke = nombre de point d'interpolation 
    {
    	
    	for (int i = 0; i < quadrature.nb_points(); i++)
    	{
    		for (int j = 0; j < quadrature.nb_points(); j++)
    		{	
    			double Somme = 0;
    			for (int pt = 0; pt < quadrature.nb_points(); pt++)
    			{
    				double eltq = 0;
    				double w = quadrature.weight(pt);
    				double k = (*coefficient)(elt_mapping.transform(quadrature.point(pt)));
    				DenseMatrix Je;
    				
    				Je = elt_mapping.jacobian_matrix(quadrature.point(pt));
    				DenseMatrix inv_Je = Je.invert_2x2();
    				DenseMatrix fin_Je = inv_Je.transpose();
    				vertex gradI = reference_functions.evaluate_grad(i, quadrature.point(pt));
    				vertex gradJ = reference_functions.evaluate_grad(j, quadrature.point(pt));
    				double detJe = elt_mapping.jacobian(quadrature.point(pt));
    				
    				vertex JexI = fin_Je.mult_2x2_2(gradI);
    				vertex JexJ = fin_Je.mult_2x2_2(gradJ);
    				
    				double dotJe = dot(JexI, JexJ);
    				
    				eltq = w * k * dotJe * detJe;
    				
    				Somme += eltq;
    			}
    			Ke.set(i,j, Somme);
    		}
    	}
    }


    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        std::cout << "Ke -> K" << '\n';
        // TODO
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (source term)" << '\n';
        // TODO
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        std::cout << "Fe -> F" << '\n';
        // TODO
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::cout << "apply dirichlet boundary conditions" << '\n';
        // TODO
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
