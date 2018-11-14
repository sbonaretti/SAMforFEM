/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VNLWRITER_H
#define VNLWRITER_H

#include <QString>

namespace vnl{

	/**
	* Abstract class for vnl writer. It sets the file name.
	*/
	class VnlWriter{
		
	public:
		
		// constructor/destructor
		VnlWriter();
		~VnlWriter();

		// accessors
		/**
		* Sets the file name
		*/
		void SetFileName (QString fileName) { _fileName = fileName;}

		// pure virtual function
		/**
		* Execute the writing
		*/
		virtual void Update() = 0;

		protected:
	
		// data members
		QString _fileName;


	};
}
#endif //VNLWRITER_H