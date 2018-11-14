/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef RENDERINGCOORDINATESYSTEM_H
#define RENDERINGCOORDINATESYSTEM_H

#include <Rendering.h>

namespace rendering{

	/**
	* It renders the origin (0,0,0) and the the versors of the three axis (x is red, y is green, z is blue). 
	* 
	* Call SetRenderer() before Update().
	*/

	class RenderingCoordinateSystem: public Rendering{

	public:

		// constructor/destructor
		RenderingCoordinateSystem();
		~RenderingCoordinateSystem();

		// overwritten virtual function
		/**
		* Executes the rendering
		*/
		void Update();

	};
}
#endif //RENDERINGCOORDINATESYSTEM_H