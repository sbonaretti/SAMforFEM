/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef RENDERINGPOINT_H
#define RENDERINGPOINT_H

#include <vtkPoints.h>
#include <Rendering.h>

namespace rendering{

	/**
	* It renders points from vtkPoints. 
	* 
	* Call SetPoints() and SetRenderer() before Update().
	* SetRadius() optional.
	*/

	class RenderingPoint: public Rendering{

	public:

		// constructor/destructor
		RenderingPoint();
		~RenderingPoint();

		// accessors
		/**
		* Sets the points to render
		*/
		void SetPoints (vtkPoints* points) { _points = points;}
		/**
		* Sets the radius of the points to render (default is 2)
		*/
		void SetRadius (double radius) { _radius = radius;}
		/**
		* Sets the interval of points to visualize (default is 1)
		*/
		void SetPointRatio (double pointRatio) { _pointRatio = pointRatio;}

				
		// overwritten virtual function
		/**
		* Executes the rendering
		*/
		void Update();

	protected:
		vtkPoints* _points;
		double _radius;
		int _pointRatio;

	};
}
#endif //RENDERINGPOINT_H