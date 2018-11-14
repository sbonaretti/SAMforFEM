/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef RENDERINGPOINTPICKER_H
#define RENDERINGPOINTPICKER_H

#include <vtkCellPicker.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>


namespace rendering{

	/**
	* It allows to pick points on the GUI. Move the mouse to the point of interest and press p.
	* It uses vtkCellPicker.
	*/
	
	class DataPacked
	{	
	public:
		DataPacked(){ red = 1.0; green = 0.0; blue = 0.0; }
		vtkRenderer* _renderer;
		vtkPolyData* _polyData;
		vtkCellPicker* _cellPicker;
		vtkPoints* _points;

		double red, green, blue;
	};

	
	class RenderingPointPicker{

	public:

		// constructor/destructor
		RenderingPointPicker();
		~RenderingPointPicker();

		// accessors
		/**
		* Sets the mesh to be picked
		*/
		void SetMesh (vtkPolyData* mesh) { _mesh = mesh;}
		/**
		* Sets the renderer
		*/
		void SetRenderer (vtkRenderer* renderer) {_renderer = renderer;}
		/**
		* Sets the color
		*/
		void SetColor (double color[3]);
		/**
		* Sets picker
		*/
		void SetPicker (vtkCellPicker* cellPicker) {_cellPicker = cellPicker;}
		/**
		* Sets double[] that will contain the picked points
		*/
		void SetPoints (vtkPoints* points) {_points = points;}
		
		// function 
		/**
		* Starts the picking
		*/
		void Update();
		/**
		* Return the picked points
		*/
		vtkPoints* GetPoints();

		
	protected:

		// accessors
		vtkPolyData* _mesh;
		vtkRenderer* _renderer;
		double _color[3];
		vtkCellPicker* _cellPicker;
		vtkPoints* _points;
	
		// data member
		DataPacked* _pack;

	};
}
#endif //RENDERINGPOINTPICKER_H