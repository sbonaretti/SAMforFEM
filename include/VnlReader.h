/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VNLREADER_H
#define VNLREADER_H

#include <QString>

namespace vnl{

	/**
	* Abstract class for vnl reader. It sets the file name.
	*/
	class VnlReader{
		
	public:
		
		// constructor/destructor
		VnlReader();
		~VnlReader();

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
#endif //VNLREADER_H