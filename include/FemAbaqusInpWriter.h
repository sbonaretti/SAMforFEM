/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef FEMABAQUSINPWRITER_H
#define FEMABAQUSINPWRITER_H

#include <vtkDoubleArray.h>
#include <vtkPoints.h> 
#include <vtkPolyData.h>

#include <QString>



namespace fem{

	/**
	* It writes the abaqus input files
	*/

	class FemAbaqusInpWriter{

	public:
		
		// constructor/destructor
		FemAbaqusInpWriter();
		~FemAbaqusInpWriter();

		// accessors
		/**
		* Sets the output file name
		*/
		void SetFileName (QString fileName) { _fileName = fileName;}
		/**
		* Sets the mesh
		*/
		void SetMesh(vtkPolyData* mesh) {_mesh = mesh;}
		/**
		* Sets the mechanical properties
		*/
		void SetMechProp(vtkDoubleArray* mechProp) {_mechProp = mechProp;}
		/**
		* Sets the boundary conditions
		*/
		void SetBoundaryConditions(vtkPoints* boundaryConditions) {_boundaryConditions = boundaryConditions;}
		/**
		* Sets the boundary condition node id
		*/
		void SetBoundaryConditionsId(vtkDoubleArray* id) {_id = id;}
		/**
		* Sets the application points
		*/
		void SetApplicationPointsID(vtkDoubleArray* appPoints) {_appPoints = appPoints;}
		
		// methods
		/**
		* 
		*/
		void WriteNodeWalkingInp();
		/**
		* 
		*/
		void WriteNodeFallingInp();
		/**
		*  for the comparison paper - fem validation
		*/
		void WriteNodeStandingInp();

	
	protected:
		
		QString _fileName;
		vtkPolyData* _mesh;
		vtkDoubleArray* _mechProp;
		vtkPoints* _boundaryConditions;
		vtkDoubleArray* _id;
		vtkDoubleArray* _appPoints;
		int _elementType;
		QString _abaqusElementType;
		double _poissonRatio; 

	};
}
#endif //FEMABAQUSINPWRITER_H