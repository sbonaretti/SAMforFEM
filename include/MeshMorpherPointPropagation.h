/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHMORPHERPOINTPROPAGATION_H
#define MESHMORPHERPOINTPROPAGATION_H

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>

#include <QStringList>

namespace mesh{
	class MeshMorpherPointPropagation{
		
	public:

		// constructor/destructor
		MeshMorpherPointPropagation();
		~MeshMorpherPointPropagation();

		// accessors
		void SetFileNames (QStringList fileNameList) {_fileNames = fileNameList;} 
		void setRenderer (vtkRenderer* renderer) {_renderer = renderer;}

		// member functions
		void Update();

	private:
	
		// data members from accessors
		QStringList _fileNames;
		vtkRenderer* _renderer;

		// functions
		void findClosestNode(double point[3], vtkPoints* points, int &minIndex);
		void writeXmlFile(QString fileName, vtkPoints* movingLandmarks, vtkDoubleArray* referenceID);

	
	};
}
#endif //MESHMORPHERPOINTPROPAGATION_H


// propagation of the 10 landmarks needed for the mesh Morpher by Ansys
	
// the 10 landmarks must be selected on the reference bone and written in a .txt file as: 
// x y z ID (each row is a landmark)
// .txt, .stl and correspondent VFs must be in the same folder

// inputs:
//   - .txt file
//   - VF (that register the images from which the meshes were created to the reference)
//   - .stl files will be automatically loaded, replacing .mhd with .stl in the file 
//     name string (they must be in the same folder)

// output:
//   - morpher .xml input file