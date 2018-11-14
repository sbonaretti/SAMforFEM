/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RenderingPointPicker.h>

#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>


namespace rendering {

	static void AddSphere(vtkObject *vtkNotUsed( caller ), unsigned long vtkNotUsed(eventId), void *sr, void *){
		
		DataPacked* data = reinterpret_cast<DataPacked*>(sr);
		
		// indeces
		int PickedCell = data->_cellPicker->GetCellId();
		if ( PickedCell == -1 )
						return;
		
		// point coords
		double coords[3];
		data->_cellPicker->GetPickPosition(coords);
		data->_points->InsertNextPoint(coords);
		
		// rendering
		vtkSphereSource* sp1=vtkSphereSource::New();
		sp1->SetRadius(2);

		vtkPolyDataMapper* sp1Mapper=vtkPolyDataMapper::New();
		sp1Mapper->SetInput(sp1->GetOutput());

		vtkActor* sp1Actor = vtkActor::New();
		sp1Actor->SetMapper(sp1Mapper);
		sp1Actor->PickableOff();
		sp1Actor->SetPosition(data->_cellPicker->GetPickPosition());
		sp1Actor->GetProperty()->SetColor(data->red, data->green, data->blue); 
			
		data->_renderer->AddActor(sp1Actor); 
 
	}


	// constructor
	RenderingPointPicker::RenderingPointPicker(){
	}

	// destructor
	RenderingPointPicker::~RenderingPointPicker(){
	}
	
	// member functions
	void RenderingPointPicker::Update(){
		    
		// DataPacked
		_pack = new DataPacked;
		_pack->_polyData = _mesh;
		_pack->_renderer = _renderer; 
		_pack->_cellPicker = _cellPicker;
		_pack->_points = _points;
		_pack->red = _color[0];
		_pack->green = _color[1];
		_pack->blue = _color[2];

		// picker
		_cellPicker->RemoveAllObservers();
			
		// setup callback for picking
		vtkCallbackCommand *cbc = vtkCallbackCommand::New();
		cbc->SetCallback(AddSphere);
		cbc->SetClientData((void *)_pack);

		_cellPicker->AddObserver( vtkCommand::EndPickEvent, cbc );
		
	}

	
	vtkPoints* RenderingPointPicker::GetPoints(){

		return _pack->_points;
	
	}

	
	void RenderingPointPicker::SetColor(double color[3]){

		_color[0] = color[0]; 
		_color[1] = color[1];
		_color[2] = color[2];
	
	}



}