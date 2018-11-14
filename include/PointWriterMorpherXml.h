/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef POINTWRITERMORPHERXML_H
#define POINTWRITERMORPHERXML_H

#include <PointWriter.h>
#include <vtkDoubleArray.h>

namespace points{

	/**
	* Writes the .xml file that is input for the ansys morpher. 
	* Use SetIdVector() to set the IDs of the reference landmarks and SetReferenceFileName() to set the reference name (that has to be written in the .xml file)
	*/
	
	class PointWriterMorpherXml: public PointWriter{
	
	public:

		// constructor/destructor
		PointWriterMorpherXml();
		~PointWriterMorpherXml();

		// member function
		/**
		* Sets the reference landmark Id vector
		*/
		void SetIdVector(vtkDoubleArray* idVector) {_idVector = idVector;}
		/**
		* Sets the reference file name (to be written in the .xml file)
		*/
		void SetReferenceFileName(QString referenceFileName) {_referenceFileName = referenceFileName;}
				
		// overwritten virtual function
		virtual void Update();


	protected:
		vtkDoubleArray* _idVector;
		QString _referenceFileName;
		
	};

}
#endif // POINTWRITERMORPHERXML_H