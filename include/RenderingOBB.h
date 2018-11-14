/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef RENDERINGOBB_H
#define RENDERINGOBB_H

#include <Rendering.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>


namespace rendering{

	/**
	* It calculates and renders the OBB of a polydata mesh (it does not render the mesh itself)
	* 
	* To calculate: call SetMesh(), ComputeOBB() before GetOBBVerteces().
	* To render: call SetMesh(), SetRenderer() and GetOBBVerteces() before Update().
	*
	* Vertex numbering: considering a vertical parallelepiped
	* - bottom base, from the front left vertex, counter-clockwise: 1-3-7-4
	* - upper base, from the front left vertex, counter-clockwise: 2-6-8-5
	
	*/

	class RenderingOBB: public Rendering{
	
	public: 

		// constructor/destructor
		RenderingOBB();
		~RenderingOBB();

		// accessors
		/**
		* Sets the mesh (not rendered)
		*/
		void SetMesh (vtkPolyData* mesh) { _mesh = mesh;}
		/**
		* Sets the points (not rendered)
		*/
		void SetPoints (vtkPoints* points) { _points = points;}
		/**
		* Gets the points (not rendered)
		*/
		vtkPoints* GetPoints() {return _points;}
		/**
		* Gets the OBB verteces
		*/
		void GetOBBverteces(double one[], double two[], double three[], double four[], 
							double five[], double six[], double seven[], double eight[] );
		/**
		* Gets the OBB center
		*/
		void GetOBBcenter(double obbCenter[]);
		/**
		* Gets the mesh center of mass
		*/
		void GetMeshCenterOfMass(double meshCenter[]);


		// member function
		
		/**
		* Computes the OBB vertices in the mesh coordinate system
		*/
		void CalculateOBB();
		/**
		* Divides the OBB in sub-OBB as needed by the class MeshBoneImplantFitting.
		* Call ComputeOBB() before ComputeOBBFree()
		*/
		void CalculateOBBtop();
		void CalculateOBBbottom();


		// overwritten virtual function
		/**
		* Executes the rendering
		*/
		void Update();

	protected:
		vtkPolyData* _mesh;
		vtkPoints* _points;
		double _one[3];
		double _two[3];
		double _three[3];
		double _four[3];
		double _five[3];
		double _six[3];
		double _seven[3];
		double _eight[3];
		double _obbCenter[3];
		double _meshCenterOfMass[3];
	
	};
}
#endif //RENDERINGOBB_H 