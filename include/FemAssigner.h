/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef FEMASSIGNER_H
#define FEMASSIGNER_H

#include <itkImage.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <QString.h>

namespace fem{

	/**
	* Abstract class for the mechanical properties assignment. 
	* 
	* In the image, voxels are taken in correspondence to the node position and converted from gray level to Young's modulus.
	* Based on Helgason et al. 2008 and Schileo et al. 2008.
	* The default law is Morgan's law (Morgan et al. 2003).
	* 
	* The assignment can be done on elements or nodes. See the classes FemAssignerElements and FemAssignerNodes.
	*/

	class FemAssigner{

	public:
		
		// image definition
		typedef float VoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Image< VoxelType, Dimension > ImageType;

		// constructor/destructor
		FemAssigner();
		~FemAssigner();

		// accessors
		/**
		* Sets the mesh file name
		*/
		void SetMesh (vtkPolyData* mesh) {_mesh = mesh;}
		/**
		* Sets the float image
		*/
		void SetImage (ImageType::Pointer image) {_image = image;}
		/**
		* Sets the first component of the assignment law
		*/	
		void SetAssignmentLawOne (double assignmentLawOne) {_assignmentLawOne = assignmentLawOne;}
		/**
		* Sets the second component of the assignment law
		*/
		void SetAssignmentLawTwo (double assignmentLawTwo) {_assignmentLawTwo = assignmentLawTwo;}
		/**
		* Gets the grey levels
		*/
		vtkDoubleArray* GetGreyLevels () {return _greyLevel;}
		/**
		* Gets the ash density
		*/
		vtkDoubleArray* GetRhoAsh () {return _rAsh;}
		/**
		* Gets the apparent density
		*/
		vtkDoubleArray* GetRhoApp () {return _rApp;}
		/**
		* Gets the Young's modulus array
		*/
		vtkDoubleArray* GetYoungModulus () {return _youngModulus;}

		// methods
		/**
		* It interpolates the grey levels surrounding the mesh nodes. 
		* The final grey level is divided by 1000 in order to convert from g/mm3 to mg/mm3
		*/
		void GreyLevelAssignmentUpdate();
		
		// pure virtual function
		/**
		* Executes the assignment
		*/
		virtual void Update() = 0;
	
		

	protected:
		vtkPolyData* _mesh;
		ImageType::Pointer _image;
		double _assignmentLawOne;
		double _assignmentLawTwo;
		vtkDoubleArray* _greyLevel;
		vtkDoubleArray* _rApp;
		vtkDoubleArray* _rAsh;
		vtkDoubleArray* _youngModulus;
	
	};
}
#endif //FEMASSIGNER_H