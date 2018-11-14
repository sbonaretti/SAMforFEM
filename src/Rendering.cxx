/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <Rendering.h>

namespace rendering {

	// constructor
	Rendering::Rendering(){
	}

	// destructor
	Rendering::~Rendering(){
	}
	
	
	// accessor
	void Rendering::SetColor(double color[3]){

		_color[0] = color[0]; 
		_color[1] = color[1];
		_color[2] = color[2];
	
	}

}