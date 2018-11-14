/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef REGISTRATIONPROCRUSTESALIGNMENT_H
#define REGISTRATIONPROCRUSTESALIGNMENT_H

#include <vtkPolyData.h>
#include <QStringList>


namespace registration{

	/**
	* It alignes meshes to their average using rigid or affine registration. 
	* It uses vtkProcrustesAlignmentFilter. 
	* Rigid registration used as default. When not specified, the similarity alignment is used.
	* (to implement: reading of differerent meshes -> flag for reading?)
	*/

	class RegistrationProcrustesAlignment{
		
	public:
		
		// constructor/destructor
		RegistrationProcrustesAlignment();
		~RegistrationProcrustesAlignment();

		// accessors
		/**
		* Sets the mesh file names
		*/	
		void SetFileNames (QStringList fileNames) { _fileNames = fileNames;}
		/**
		* Sets the registration type. 0 is rigid alignment, 1 is affine alignment.
		*/
		void SetRegistrationFlag (int flag) {_flag = flag;}
		/**
		* Sets the source mesh
		*/	
		void SetSourceMesh (vtkPolyData* sourceMesh) { _sourceMesh = sourceMesh;}
		/**
		* Sets the target mesh
		*/	
		void SetTargetMesh (vtkPolyData* targetMesh) { _targetMesh = targetMesh;}
		/**
		* Get the source mesh alligned
		*/	
		vtkPolyData* GetAlignedMesh() {return _outputMesh;}



		// member functions
		/**
		* Execute the alignment
		*/	
		void Update();
		/**
		* Execute the alignment using ICP - useful for meshes with different number of points - only for 2 meshes
		* it sets 2 vtkPointsSets
		*/	
		void IterativeClosestPoints();
	
	protected:
	
		// data members from accessors
		QStringList _fileNames;
		int _flag;
		// for ICP
		vtkPolyData* _sourceMesh;
		vtkPolyData* _targetMesh;
		vtkPolyData* _outputMesh;
	
	};
}
#endif //REGISTRATIONPROCRUSTESALIGNMENT_H