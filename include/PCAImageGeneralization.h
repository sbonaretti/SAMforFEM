/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef PCAIMAGEGENERALIZATION_H
#define PCAIMAGEGENERALIZATION_H

#include <QStringList>

#include <PCAImageValidation.h>

namespace pca{

	/**
	* It calculates the model ability to be general.
	* G = 1/N * sum on i (1->N) min ||instance a - dataset i|| (N = dataset dimension).
	* By Davies R., Statistical model of shape/ pag. 238.
	*/
	
	class PCAImageGeneralization: public PCAImageValidation{
	
	
	public:
		// constructor/destructor
		PCAImageGeneralization();
		~PCAImageGeneralization();

		// member function
		/**
		* Sets the file names of the created dataset
		*/
		void SetInstanceFileNames (QStringList instanceFileNames) {_instanceFileNames = instanceFileNames;}

		// overwritten virtual function
		/**
		* Computes the calculation
		*/
		virtual void Update();

	protected:
		QStringList _instanceFileNames;

	};

}
#endif // PCAIMAGEGENERALIZATION_H

