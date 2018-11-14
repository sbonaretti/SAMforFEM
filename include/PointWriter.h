/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef POINTWRITER_H
#define POINTWRITER_H

#include <vtkPoints.h>
#include <QString>

#include <QProcess.h>
#include <QFile.h>
#include <QTextStream.h>


namespace points{

	/**
	* Abstract class for point writer. It sets the file name and the points.
	*/

	class PointWriter{
	
	public:

		// constructor/destructor
		PointWriter();
		~PointWriter();

		// accessors
		/**
		* Sets the file name
		*/
		void SetFileName (QString fileName) { _fileName = fileName;}
		/**
		* Sets the points
		*/
		void SetInput (vtkPoints* points) {_points = points;}
				
		// pure virtual function
		/**
		* Execute the writing
		*/
		virtual void Update() = 0;
		

	protected:

		// data members from accessors
		QString _fileName;
		vtkPoints* _points;

	};

}
#endif // POINTWRITER_H