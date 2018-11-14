/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHQUALITY_H
#define MESHQUALITY_H

#include <vtkPolyData.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <QStringList>


namespace mesh{

	/**
	* Evaluation of mesh quality.
	* Mesh considered are tetrahedron ones, both linear and quadratic (treated as linear - only vertices considered)
	* It uses vtkMeshQuality which uses TetMetric.cpp which is part of VERDICT (vtk/Utilities/Verdict)
	* It write files containing the final results of the calculations
	* 
	* CalculateJacobian() calculates the scaled Jacobian of each element ( minimum of the jacobian divided by the lengths of 3 edge vectors)
	* CalculateEdgeRatio() calculates the edge ratio of each element (longest edge / shortest edge)
	* CalculateMinAngle() calculates the min angle of each element
	*/

	class MeshQuality{
		
	public:
		
		// constructor/destructor
		MeshQuality();
		~MeshQuality();

		// accessors
		/**
		* Sets the input mesh filenames
		*/	
		void SetMeshFileNames (QStringList meshFileNames) { _meshFileNames = meshFileNames;}
		/**
		* Gets the jacobian of each tetra
		*/	
		vnl_vector<double> GetJacobian() {return _jacobians;} 
		/**
		* Gets the edge ratio of each tetra
		*/	
		vnl_vector<double> GetEdgeRatio() {return _edgeRatio;} 
		/**
		* Gets the minumum angle of each tetra
		*/	
		vnl_vector<double> GetMinAngle() {return _minAngle;} 

		
		// member functions
		/**
		* Calculates the Jacobian of each element
		*/	
		void CalculateJacobian();
		/**
		* Calculates the edge ratio of each element
		*/	
		void CalculateEdgeRatio();
		/**
		* Calculates the minimum angle of each element
		*/	
		void CalculateMinAngle();


	protected:
		// data members from accessors
		QStringList _meshFileNames;
		/**
		* Contains the jacobians of all elements for all meshes
		*/	
		vnl_vector<double> _jacobians;
		/**
		* Contains the number of null or negative jacobians
		*/
		vnl_matrix<double> _jacobiansQuality;
		/**
		* Contains the edge ratios of all elements for all meshes
		*/	
		vnl_vector<double> _edgeRatio;
		/**
		* Contains mesh id, edge ratio average and edge ratio standard deviation 
		*/	
		vnl_matrix<double> _edgeRatioQuality;
		/**
		* Contains the min angles of all elements for all meshes
		*/	
		vnl_vector<double> _minAngle;
		/**
		* Contains mesh id, min angle average and min angle standard deviation 
		*/	
		vnl_matrix<double> _minAngleQuality;


	};
}
#endif //MESHQUALITY_H