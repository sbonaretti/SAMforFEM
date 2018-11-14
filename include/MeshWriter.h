/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHWRITER_H
#define MESHWRITER_H

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>

#include <QString>


namespace mesh{

	/**
	* Abstract class for mesh writer. It sets the file name and the mesh.
	*/
	class MeshWriter{
		
	public:
		
		// constructor/destructor
		MeshWriter();
		~MeshWriter();

		// accessors
		/**
		* Sets the file name
		*/
		void SetFileName (QString fileName) { _fileName = fileName;}
		/**
		* Sets the mesh flag on in order to write the mesh in the file
		*/
		void MeshOn(){_meshFlag = true;}
		/**
		* Sets the mesh
		*/
		void SetMesh (vtkPolyData* mesh) {_mesh = mesh;}
		/**
		* Sets the Young's modulus flag on in order to write the elements mechanical properties
		*/
		void YoungsModulusElemOn(){_youngModulusElemFlag = true;}
		/**
		* Sets the Young's modulus flag on in order to write the node mechanical properties
		*/
		void YoungsModulusNodeOn(){_youngModulusNodeFlag = true;}
		/**
		* Sets Young's modulus
		*/
		void SetYoungsModulus (vtkDoubleArray* youngModulus) {_youngModulus = youngModulus;}
		/**
		* Setst the load flag on
		*/
		void SetLoadOn() {_loadFlag = true;}
		/**
		* Sets the load magnitude array
		*/
		void SetLoadMagnitude(vtkDoubleArray* loadMagnitude) {_loadMagnitude = loadMagnitude;}
		/**
		* Sets the load application points
		*/
		void SetLoadPoints(vtkPoints* loadPoints) {_loadPoints = loadPoints;}
		/**
		* Sets the boundary condition flag on
		*/
		void SetBoundaryConditionOn(){_bcFlag = true;}
		/**
		* Sets the boundary condition types (1=encastre)
		*/
		void SetBoundaryConditionType(vtkDoubleArray* bcType){_bcType = bcType;}
		/**
		* Sets the boundary condition application points
		*/
		void SetBoundaryConditionPoints(vtkPoints* bcPoints){_bcPoints = bcPoints;}
		


		// pure virtual function
		/**
		* Execute the writing
		*/
		virtual void Update() = 0;
	
	
	protected:
	
		// data members
		QString _fileName;
		bool _meshFlag;
		vtkPolyData* _mesh;
		bool _youngModulusElemFlag;
		bool _youngModulusNodeFlag;
		vtkDoubleArray* _youngModulus;
		bool _loadFlag;
		vtkDoubleArray* _loadMagnitude;
		vtkPoints* _loadPoints;
		bool _bcFlag;
		vtkDoubleArray* _bcType;
		vtkPoints* _bcPoints;


		

	};
}
#endif //MESHWRITER_H