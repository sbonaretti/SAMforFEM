/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MESHREADER_H
#define MESHREADER_H

#include <vtkPolyData.h>
#include <QString>


namespace mesh{

	/**
	* Abstract class for mesh reader. It sets the file name and gets the mesh as output.
	*/
	class MeshReader{
		
	public:
		
		// constructor/destructor
		MeshReader();
		~MeshReader();

		// accessors
		/**
		* Sets the file name
		*/
		void SetFileName (QString fileName) { _fileName = fileName;}
		/**
		* Gets the mesh
		*/
		vtkPolyData* GetOutput () {return _mesh;}
		
		// pure virtual function
		/**
		* Executes the reading
		*/
		virtual void Update() = 0;
		
	
	protected:
	
		// data members from accessors
		QString _fileName;
		vtkPolyData* _mesh;

	};
}
#endif //MESHREADER_H