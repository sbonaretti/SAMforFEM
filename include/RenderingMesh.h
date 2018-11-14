/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef RENDERINGMESH_H
#define RENDERINGMESH_H

#include <Rendering.h>

#include <vtkPolyData.h>


namespace rendering{

	/**
	* It renders meshes from vtkPolyData. For tetra quadratic meshes use RenderingQuadratiTetraMesh.
	* 
	* Call SetMesh() and SetRenderer() before Update().
	*/

	class RenderingMesh: public Rendering{
	
	public: 

		// constructor/destructor
		RenderingMesh();
		~RenderingMesh();

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
#endif //RENDERINGMESH_H 