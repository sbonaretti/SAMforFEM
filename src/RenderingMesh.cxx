/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RenderingMesh.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>


namespace rendering{
	
	// constructor
	RenderingMesh::RenderingMesh(){
		
		_color[0] = 232.0/255.0; //default color
		_color[1] = 184.0/255.0;
		_color[2] = 45.0/255.0;

	}

	
	// destructor
	RenderingMesh::~RenderingMesh(){
	}
	
	
	// overwritten virtual function
	void RenderingMesh::Update(){
					
		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInput(_mesh);
		mapper->ScalarVisibilityOff();
		
		vtkActor* actor = vtkActor::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(_color);
		
		_renderer->AddActor(actor);
		_renderer->ResetCamera();
		_renderer->GetRenderWindow()->Render();
	
	}

}