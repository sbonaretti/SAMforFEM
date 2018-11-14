/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RenderingQuadraticTetraMesh.h>

#include <vtkActor.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkUnstructuredGrid.h>


namespace rendering{
	
	// constructor
	RenderingQuadraticTetraMesh::RenderingQuadraticTetraMesh(){
		
		_color[0] = 232.0/255.0; //default color
		_color[1] = 184.0/255.0;
		_color[2] = 45.0/255.0;

	}

	
	// destructor
	RenderingQuadraticTetraMesh::~RenderingQuadraticTetraMesh(){
	}
	
	
	// overwritten virtual function
	void RenderingQuadraticTetraMesh::Update(){

		vtkUnstructuredGrid* mesh = vtkUnstructuredGrid::New();
		mesh->SetPoints(_mesh->GetPoints());
		mesh->SetCells(VTK_QUADRATIC_TETRA, _mesh->GetPolys());

		vtkGeometryFilter* extractOuterSurface = vtkGeometryFilter::New();
		extractOuterSurface->SetInput(mesh);
		extractOuterSurface->MergingOff();
		extractOuterSurface->Update();

		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInput(extractOuterSurface->GetOutput());
		mapper->ScalarVisibilityOff();
		
		vtkActor* actor = vtkActor::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(_color);
		
		_renderer->AddActor(actor);
		_renderer->ResetCamera();
		_renderer->GetRenderWindow()->Render();


		}
}