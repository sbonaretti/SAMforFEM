/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef FEMFORCE_H
#define FEMFORCE_H

#include <vtkDoubleArray.h>
#include <vtkPoints.h> 
#include <vtkPolyData.h> 

#include <itkImage.h> 


namespace fem{

	/**
	* 
	*/

	class FemForce{

	public:
		
		// constructor/destructor
		FemForce();
		~FemForce();

		// 3D image
		typedef signed short VoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Image< VoxelType, Dimension > ImageType;
		
		// vector image
		typedef float FieldVoxelType;
		typedef itk::Vector<FieldVoxelType,Dimension> VectorType;
		typedef itk::Image<VectorType,Dimension> FieldType; 


		// accessors
		/**
		* Sets the application points
		*/
		void SetPoints(vtkPoints* points) {_points = points;}
		/**
		* Sets the application points id
		*/
		void SetPointsId(vtkDoubleArray* idVector) {_idVector = idVector;}
		/**
		* Sets the mesh
		*/
		void SetMesh(vtkPolyData* mesh) {_mesh = mesh;}
		/**
		* Sets the stationary velocity field
		*/
		void SetSVF(FieldType::Pointer field) {_field = field;}		
		/**
		* Sets the body weight
		*/
		void SetBodyWeight(double bodyWeight) {_bodyWeight = bodyWeight;}	
		
		/**
		* Gets the application or propagated points
		*/
		vtkPoints* GetPoints() {return _points;}
		/**
		* Gets the mesh ID of the propagated points (ID starts from 1, not 0)
		*/
		vtkDoubleArray* GetPointsID() {return _idVector;}

		/**
		* Gets the force (decomposed in the three directions) 
		*/
		double GetForce() {return _force[3];}
		
		

		// methods
		/**
		* It calculates the forces that represent walking (Heller et al. 2005)
		* Call SetPoints()
		* The input point array contains: 
		* - P0: above fovea 
		* - P1: distal posterior greater trochanter (bottom of the bump)
		* - P2: distal posterior greater trochanter
		* - F1: proximal middle greater trochanter (middle of the crest)
		* - F3: distal intercondyle
		* - F5: lateral condyle (middle of the bump)
		*
		* The output point array contains:
		* - hip contact in P0 (three load components)
		* - abductor force in P1 (three load components)
		* - proximal tensor fascia latae force in P1 (three load components)
		* - distal tensor fascia latae force in P1 (three load components)
		* - vastus lateralis in P2 (three load components)
		* - point A as versor of F3-P0 (for abaqus coordinate system - boundary conditions)
		* - point B as versor of F5-P0 (for abaqus coordinate system - boundary conditions)
		* - point C as versor of F3-F1 (for abaqus coordinate system - forces)
		* - point D as versor of F5-F1 (for abaqus coordinate system - forces)
		* - femur length, body height, BMI
		* The output id array contains:
		* - node id of P0
		* - node id of P1
		* - node id of P1
		* - node id of P1
		* - node id of P2
		* - 0 (fictitious value) A 
		* - 0 (fictitious value) B
		* - 0 (fictitious value) C
		* - 0 (fictitious value) D
		* - body weight
		* (if written with PointWriterXyzId.h the array becomes the fourth column of the file)
		*/
		void WalkingFromFemulLenght();
		/**
		* Recalls the body weight and define the forces for the falling situation
		*/
		void SideFalling();
		/**
		* Force directed from P0 to F3. Same input point array as WalkingFromFemulLenght()
		*/
		void Standing();
				
		/**
		* Propagates the force and boundary condition application points
		* Call SetPoints(), SetMesh(), SetSVF()
		*/
		void PropagateApplicationPoints();
		
	
		/**
		* Find the mesh closest node to the given point 
		* minIndex is the closest point index; closestX, closestY, closestZ are the x,y,z closest point coordinates 
		*/
		void findClosestNode(vtkPolyData* mesh, double point[3], int &minIndex, double &closestX, double &closestY, double &closestZ);
				

	protected:

		
		vtkPoints* _points;
		vtkPolyData* _mesh;
		FieldType::Pointer _field;
		vtkDoubleArray* _idVector;
		double _force[3];
		double _bodyWeight;
	};
}
#endif //FEMFORCE_H