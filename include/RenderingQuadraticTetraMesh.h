/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef RENDERINGQUADRATICTETRAMESH_H
#define RENDERINGQUADRATICTETRAMESH_H

#include <Rendering.h>

#include <vtkPolyData.h>


namespace rendering{

	/**
	* It renders a quadratic tetra mesh. It takes the volume mesh as polydata, tranforms it to
	* vtkUnstrucuredGrid, applies the vtkGeometryFilter that works for rendering (not to extract the surface mesh).
	* For linear or quadratic triangular surface meshes use RenderingMesh.
	* 
	* Call SetMesh() and SetRenderer() before Update().
	*/

	class RenderingQuadraticTetraMesh: public Rendering{
	
	public: 

		// constructor/destructor
		RenderingQuadraticTetraMesh();
		~RenderingQuadraticTetraMesh();

		// accessors
		/**
		* Sets the mesh to render
		*/
		void SetMesh (vtkPolyData* mesh) { _mesh = mesh;}
				
		// overwritten virtual function
		/**
		* Executes the rendering
		*/
		void Update();

	protected:
		vtkPolyData* _mesh;
		
	};
}
#endif //RENDERINGQUADRATICTETRAMESH_H 