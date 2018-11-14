/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef RENDERINGPOINTWITHCOLORBAR_H
#define RENDERINGPOINTWITHCOLORBAR_H

#include <vtkDoubleArray.h>
#include <vtkPoints.h>

#include <RenderingPoint.h>

namespace rendering{

	/**
	* It renders points from vtkPoints with a colorbar(vtkDoubleArray). 
	* 
	* Call SetPoints() SetColorbar() and SetRenderer() before Update().
	* SetRadius() optional. 
	* SetColorRange() sets the min and max colors of the color bar. Default is b/w.
	*/

	class RenderingPointWithColorbar: public RenderingPoint{

	public:

		// constructor/destructor
		RenderingPointWithColorbar();
		~RenderingPointWithColorbar();

		// accessors
		/**
		* Sets the colorbar
		*/
		void SetColorBar (vtkDoubleArray* colorBar) { _colorBar = colorBar;}
		/**
		* Sets the colorbar range
		*/
		void SetColorRange (double colorRange[2]); 

		
		// overwritten virtual function
		/**
		* Executes the rendering (b/w)
		*/
		void UpdateBW();
		/**
		* Executes the rendering (colors)
		*/
		void UpdateC();
		/**
		* Executes the rendering, labelling the bar as E(GPa)
		*/
		void UpdateYoungModulus();

	protected:
		vtkDoubleArray* _colorBar;
		double _colorRange[2];

	};
}
#endif //RENDERINGPOINTWITHCOLORBAR_H