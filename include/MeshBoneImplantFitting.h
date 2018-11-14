/*
 * ISTB - University of Bern, Mauricio Reyes, Serena Bonaretti
 */

#ifndef MESHBONEIMPLANTFITTING_H
#define MESHBONEIMPLANTFITTING_H

#include <vtkPolyData.h>

namespace mesh{

	/**
	* 
	*/
	class MeshBoneImplantFitting{
		
	public:
		
		// constructor/destructor
		MeshBoneImplantFitting();
		~MeshBoneImplantFitting();

		// accessors
		/**
		* Sets the bone mesh
		*/
		void SetBoneMesh(vtkPolyData* boneMesh){_boneMesh = boneMesh;}
		/**
		* Sets the implant mesh
		*/
		void SetImplantMesh(vtkPolyData* implantMesh){_implantMesh = implantMesh;}
		/**
		* Gets the bone mesh
		*/
		vtkPolyData* GetBoneMesh() {return _boneMesh;}		
		/**
		* Gets the implant mesh
		*/
		vtkPolyData* GetImplantMesh() {return _implantMesh;}		
		/**
		* Sets the y resolution
		*/
		void SetYResolution(double yResolution) {_yResolution = yResolution;}
		/**
		* Sets the z resolution
		*/
		void SetZResolution(double zResolution) {_zResolution = zResolution;}
		/**
		* Sets the y step
		*/
		void SetYStep(double yStep) {_yStep = yStep;}
		/**
		* Sets the z step
		*/
		void SetZStep(double zStep) {_zStep = zStep;}
		/**
		* 
		*/
		void SetImplantPosition(int positionFlag) {_positionFlag = positionFlag;}
		



		// member functionds
		/**
		* Executes the fitting
		*/
		void Update();
		/**
		* 
		*/
		void ImplantPositionUpdate();
		
	protected:
	
		// data members
		vtkPolyData* _boneMesh;
		vtkPolyData* _implantMesh;
		double _yResolution;
		double _zResolution;
		double _yStep;
		double _zStep;
		double _positionFlag;

	private:
		/**
		* It translates the given mesh of the amount dx, dy, dz and returns the mesh back 
		*/
		void TranslateMesh(vtkPolyData* mesh, double dx,double dy,double dz);
		/**
		* 
		*/
		void RotateMesh(vtkPolyData* mesh, double oneI[], double twoI[], double threeI[], double fourI[],
										   double oneB[], double twoB[], double threeB[], double fourB[]);
		/**
		* 
		*/
		//void CenterOfMass(vtkPolyData *data,double com[3]);
		/**
		* 
		*/
		//void CheckCollision(vtkPolyData *implantPoly,vtkPolyData *bonePoly, vtkPoints* distances);
		/**
		* 
		*/
		//void InVolume(vtkPolyData *pd,double point[3], double& x, double& y, double& z);
		/**
		* 
		*/
		//void OBB(vtkPolyData* mesh, )	

	};
}
#endif //MESHWRITER_H