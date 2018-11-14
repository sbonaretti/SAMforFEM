/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef PCAIMAGESPECIFICITY_H
#define PCAIMAGESPECIFICITY_H

#include <QStringList>

#include <PCAImageValidation.h>

namespace pca{
	class PCAImageSpecificity: public PCAImageValidation{
		
	public:
		// constructor/destructor
		PCAImageSpecificity();
		~PCAImageSpecificity();

		// member function
		void SetInstanceFileNames (QStringList instanceFileNames) {_instanceFileNames = instanceFileNames;}

		// overwritten virtual function
		virtual void Update();

	protected:
		QStringList _instanceFileNames;

	};

}
#endif PCAIMAGESPECIFICITY_H