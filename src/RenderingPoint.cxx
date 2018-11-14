/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RenderingPoint.h>

#include <vtkActor.h>
#include <vtkGlyph3D.h>
#include <vtkMaskPoints.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkSphereSource.h>


namespace rendering{

	// constructor
	RenderingPoint::RenderingPoint(){
			
		_color[0] = 1.0; //default color
		_color[1] = 0.0;
		_color[2] = 0.0;

		_radius = 2;
		_pointRatio = 1;
		
	}

	// destructor
	RenderingPoint::~RenderingPoint(){
	}

	// overwritten virtual function
	void RenderingPoint::Update(){

		vtkPolyData* points = vtkPolyData::New();
		points->SetPoints(_points);
		
		vtkMaskPoints* maskPoints = vtkMaskPoints::New();
		maskPoints->SetInput(points);
		maskPoints->SetOnRatio(_pointRatio); // interval of points to visualize (one every...)
		
		vtkSphereSource* sphere = vtkSphereSource::New();
		sphere->SetRadius(_radius); // sphere diameter
		
		vtkGlyph3D* glyph = vtkGlyph3D::New();
		glyph->SetInput(maskPoints->GetOutput());
		glyph->SetSource(sphere->GetOutput());
		
		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInput(glyph->GetOutput());
		
		vtkActor* actor = vtkActor::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(_color);
		
		_renderer->AddActor(actor);
		_renderer->ResetCamera();
		_renderer->GetRenderWindow()->Render();



	}
}