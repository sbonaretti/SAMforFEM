/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef POINTREADER_H
#define POINTREADER_H

#include <vtkPoints.h>
#include <QString>


namespace points{
	class PointReader{
	
	public:

		// constructor/destructor
		PointReader();
		~PointReader();

		// accessors
		void SetFileName (QString fileName) { _fileName = fileName;}
		vtkPoints* GetOutput() { return _points;}
				
		// pure virtual function
		virtual void Update() = 0;
		

	protected:

		// data members from accessors
		QString _fileName;
		vtkPoints* _points;

	};

}
#endif // POINTREADER_H