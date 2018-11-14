/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef RENDERING_H
#define RENDERING_H

#include <vtkRenderer.h>
#include <QString>


namespace rendering{

	/**
	* Abstract class for rendering. It sets the color of the object to render and the renderer.
	*/

	class Rendering{

	public:

		// constructor/destructor
		Rendering();
		~Rendering();

		// accessors
		/**
		* Sets the color
		*/
		void SetColor(double color[3]);
		/**
		* Sets the renderer
		*/
		void SetRenderer (vtkRenderer* renderer) {_renderer = renderer;}

		// pure virtual functions
		/**
		* Executes the rendering
		*/
		virtual void Update()=0;

	protected:
		QString _fileName;
		double _color[3];
		vtkRenderer* _renderer;
	
	};
}
#endif //RENDERING_H