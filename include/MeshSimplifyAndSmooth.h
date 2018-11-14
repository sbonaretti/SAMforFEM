/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHSIMPLIFYANDSMOOTH_H
#define MESHSIMPLIFYANDSMOOTH_H

#include <vtkPolyData.h>
#include <QString>


namespace mesh{


	/**
	* Simplification and smoothing of a surface mesh.
	* Simplification reduces the number of points. Use SetSimplifyValue() to define the percentage of reduction.
	* Smoothing smooths the mesh. Use SetSmoothIteration () to define the number of iterations.
	* ! After the use of this class, the mesh can have lost its topology. MRFSurface is used to rebuild it (path specified in Update())
	*/
	class MeshSimplifyAndSmooth{
		
	public:
		
		// constructor/destructor
		MeshSimplifyAndSmooth();
		~MeshSimplifyAndSmooth();

		// accessors
		/**
		* Sets the input mesh
		*/	
		void SetInput (vtkPolyData* mesh) { _mesh = mesh;}
		/**
		* Sets the simplification percentage. Use a number between 0.0 and 1.0
		*/
		void SetSimplifyValue (double value) { _simplifyValue = value;}
		/**
		* Sets the smoothing iterations number
		*/
		void SetSmoothIteration (int iteration) { _smoothNofIteration = iteration;}
		/**
		* Gets the output mesh
		*/
		vtkPolyData* GetOutput () {return _mesh;}
		

		// member functions
		/**
		* Execute the simplify and smoothing
		*/	
		void Update();
		
	
	protected:
	
		// data members from accessors
		vtkPolyData* _mesh;
		double _simplifyValue;
		int _smoothNofIteration;

	};
}
#endif //MESHSIMPLIFYANDSMOOTH_H