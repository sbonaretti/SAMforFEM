/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RenderingCoordinateSystem.h>

#include <vtkActor.h>
#include <vtkLineSource.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkSphereSource.h>


namespace rendering{

	// constructor
	RenderingCoordinateSystem::RenderingCoordinateSystem(){
	}

	// destructor
	RenderingCoordinateSystem::~RenderingCoordinateSystem(){
	}

	// overwritten virtual function
	void RenderingCoordinateSystem::Update(){
		
		// variables
		double or[3]; or[0] = 0.0; or[1] = 0.0; or[2] = 0.0;
		double radius = 5.0;
		double ptX[3]; ptX[0] = 100.0; ptX[1] = 0.0; ptX[2] = 0.0;  
		double ptY[3]; ptY[0] = 0.0; ptY[1] = 100.0; ptY[2] = 0.0;  
		double ptZ[3]; ptZ[0] = 0.0; ptZ[1] = 0.0; ptZ[2] = 100.0;  
		double lineWidth = 5.0;

		// origin
		vtkSphereSource* origin = vtkSphereSource::New();
		origin->SetCenter(or);
		origin->SetRadius(radius);
		vtkPolyDataMapper* mapperC = vtkPolyDataMapper::New();
		mapperC->SetInput(origin->GetOutput());
		vtkActor* actorC = vtkActor::New();
		actorC->SetMapper(mapperC);
		_color[0] = 0.0; _color[1] = 0.0; _color[2] = 0.0;
		actorC->GetProperty()->SetColor(_color);
		_renderer->AddActor(actorC);
		
		// axis X
		vtkLineSource* lineX = vtkLineSource::New();
		lineX->SetPoint1(or);
		lineX->SetPoint2(ptX);
		vtkPolyDataMapper * mapperX = vtkPolyDataMapper ::New();
		mapperX->SetInput(lineX->GetOutput());
		vtkActor * actorX = vtkActor ::New();
		actorX->SetMapper(mapperX);
		_color[0] = 1.0; _color[1] = 0.0; _color[2] = 0.0;
		actorX->GetProperty()->SetColor(_color);
		actorX->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actorX);
		
		// axis Y
		vtkLineSource* lineY = vtkLineSource::New();
		lineY->SetPoint1(or);
		lineY->SetPoint2(ptY);
		vtkPolyDataMapper* mapperY = vtkPolyDataMapper::New();
		mapperY->SetInput(lineY->GetOutput());
		vtkActor* actorY = vtkActor::New();
		actorY->SetMapper(mapperY);
		_color[0] = 0.0; _color[1] = 1.0; _color[2] = 0.0;
		actorY->GetProperty()->SetColor(_color);
		actorY->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actorY);
		
		// axis Z
		vtkLineSource* lineZ = vtkLineSource::New();
		lineZ->SetPoint1(or);
		lineZ->SetPoint2(ptZ);
		vtkPolyDataMapper* mapperZ = vtkPolyDataMapper::New();
		mapperZ->SetInput(lineZ->GetOutput());
		vtkActor* actorZ = vtkActor::New();
		actorZ->SetMapper(mapperZ);
		_color[0] = 0.0; _color[1] = 0.0; _color[2] = 1.0;
		actorZ->GetProperty()->SetColor(_color);
		actorZ->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actorZ);
		
		_renderer->ResetCamera();
		_renderer->GetRenderWindow()->Render();

		// cleaning up
		origin->Delete();
		mapperC->Delete();
		actorC->Delete();
		lineX->Delete();
		mapperX->Delete();
		actorX->Delete();
		lineY->Delete();
		mapperY->Delete();
		actorY->Delete();
		lineZ->Delete();
		mapperZ->Delete();
		actorZ->Delete();

	}


}

