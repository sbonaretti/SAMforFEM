/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RenderingPointWithColorbar.h>

#include <vtkActor.h>
#include <vtkCoordinate.h>
#include <vtkGlyph3D.h>
#include <vtkLookupTable.h>
#include <vtkMaskPoints.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkScalarBarActor.h>
#include <vtkSphereSource.h>
#include <vtkTextProperty.h>

#include <vtkIntArray.h>
#include <vtkColorTransferFunction.h>
#include <vtkPointData.h>



namespace rendering{

	// constructor
	RenderingPointWithColorbar::RenderingPointWithColorbar(){

		_colorRange[0] = 0.0; 
		_colorRange[1] = 1.0;
	}

	// destructor
	RenderingPointWithColorbar::~RenderingPointWithColorbar(){
	}

	// accessor
	void RenderingPointWithColorbar::SetColorRange(double colorRange[2]){

		_colorRange[0] = colorRange[0]; 
		_colorRange[1] = colorRange[1];
	
	}

	// overwritten virtual function
	void RenderingPointWithColorbar::UpdateBW(){

		vtkPolyData* polyData = vtkPolyData::New();
		polyData->SetPoints(_points);

		double range[2]; // range for visualization
		_colorBar->GetRange(range);

		// add colors to vtkPolyData
		vtkIntArray* colors = vtkIntArray::New();
		colors->DeepCopy(_colorBar);
		colors->SetName("col");
		polyData->GetPointData()->AddArray(colors);

		
		vtkMaskPoints* maskPoints = vtkMaskPoints::New();
		maskPoints->SetInput(polyData);
		maskPoints->SetOnRatio(_pointRatio); // interval of points to visualize (one every...)
		
		vtkSphereSource* sphere = vtkSphereSource::New();
		sphere->SetRadius(_radius); // sphere diameter
		sphere->SetPhiResolution(4); // sphere resolution
		sphere->SetThetaResolution(4);

		vtkGlyph3D* glyph = vtkGlyph3D::New();
		glyph->SetInput(maskPoints->GetOutput());
		glyph->SetSource(sphere->GetOutput());

		// lookup table (converts numbers to colors)
		vtkLookupTable* lut = vtkLookupTable::New();
		if (_colorRange[0] == 0.0 && _colorRange[1] == 1){
			lut->SetTableRange (range[0], range[1]);
			lut->SetValueRange (_colorRange[0], _colorRange[1]);
			lut->SetSaturationRange(0.0, 0.0); // no color saturation
			lut->SetRampToLinear();
			lut->SetNumberOfColors(256);
			lut->Build();
		}
		else {
			lut->SetTableRange (range[0], range[1]);
			lut->SetHueRange (_colorRange[0], _colorRange[1]);
			lut->Build();
			lut->SetNumberOfColors(256);
		}
		
		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInput(glyph->GetOutput());
		mapper->SetLookupTable(lut);
		mapper->SetScalarRange(lut->GetTableRange());
		mapper->ScalarVisibilityOn();
		mapper->SetScalarModeToUsePointFieldData(); // SetScalarModeToUsePointData() doesnt allow the points to be colored!!!
		mapper->ColorByArrayComponent("col",0);
		mapper->SetColorModeToMapScalars();
		
		vtkActor* actor = vtkActor::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(_color);

		// colorbar
		vtkScalarBarActor* scalarBar = vtkScalarBarActor::New();
		//scalarBar->SetTitle("test");
		scalarBar->SetNumberOfLabels(5);
		//scalarBar->SetOrientationToHorizontal();
		scalarBar->SetOrientationToVertical();
		// label text
		scalarBar->GetLabelTextProperty()->SetColor(51.0/255.0, 51.0/255.0, 51.0/255.0);
		scalarBar->GetLabelTextProperty()->ShadowOff();
		scalarBar->GetLabelTextProperty()->SetFontSize(25);
		// title text
		scalarBar->GetTitleTextProperty()->SetColor(51.0/255.0, 51.0/255.0, 51.0/255.0);
		scalarBar->GetTitleTextProperty()->ShadowOff();
		scalarBar->GetLabelTextProperty()->SetFontSize(25);
		scalarBar->SetLookupTable(mapper->GetLookupTable());
		
		vtkCoordinate* barCoord = vtkCoordinate::New(); // colorbar position in window and dimensions
		barCoord = scalarBar->GetPositionCoordinate();
		barCoord->SetCoordinateSystemToNormalizedViewport();
		barCoord->SetValue(0.1,0.05);
		//scalarBar->SetWidth(0.8); // when the bar is horizontal
		//scalarBar->SetHeight(0.1);
		scalarBar->SetWidth(0.1); // when the bar is vertical
		scalarBar->SetHeight(1.0);
		scalarBar->SetLookupTable(lut);
		
		_renderer->AddActor(actor);
		_renderer->AddActor2D(scalarBar);
		_renderer->ResetCamera();
		_renderer->GetRenderWindow()->Render();

		// cleaning up
		/*
		polyData->Delete();
		colors->Delete();
		maskPoints->Delete();
		sphere->Delete();
		glyph->Delete();
		lut->Delete();
		mapper->Delete();
		actor->Delete();
		scalarBar->Delete();
		barCoord->Delete();
		*/
	
	}

	// overwritten virtual function
	void RenderingPointWithColorbar::UpdateC(){

		vtkPolyData* polyData = vtkPolyData::New();
		polyData->SetPoints(_points);

		double range[2]; // range for visualization
		_colorBar->GetRange(range);

		// add colors to vtkPolyData
		vtkIntArray* colors = vtkIntArray::New();
		colors->DeepCopy(_colorBar);
		colors->SetName("col");
		polyData->GetPointData()->AddArray(colors);

		
		vtkMaskPoints* maskPoints = vtkMaskPoints::New();
		maskPoints->SetInput(polyData);
		maskPoints->SetOnRatio(5); // interval of points to visualize (one every...)
		
		vtkSphereSource* sphere = vtkSphereSource::New();
		sphere->SetRadius(0.5); // sphere diameter
		
		vtkGlyph3D* glyph = vtkGlyph3D::New();
		glyph->SetInput(maskPoints->GetOutput());
		glyph->SetSource(sphere->GetOutput());

		// lookup table (converts numbers to colors)
		vtkLookupTable* lut = vtkLookupTable::New();
		lut->SetTableRange (range[0], range[1]);
		//lut->SetHueRange (0.0, 0.67); // red = min; blue = max
		lut->SetHueRange (0.67, 0.0); // blue = min; red = max
		lut->Build();
		lut->SetNumberOfColors(256);

		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInput(glyph->GetOutput());
		mapper->SetLookupTable(lut);
		mapper->SetScalarRange(lut->GetTableRange());
		mapper->ScalarVisibilityOn();
		mapper->SetScalarModeToUsePointFieldData(); // SetScalarModeToUsePointData() doesnt allow the points to be colored!!!
		mapper->ColorByArrayComponent("col",0);
		mapper->SetColorModeToMapScalars();
		
		vtkActor* actor = vtkActor::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(_color);

		// colorbar
		vtkScalarBarActor* scalarBar = vtkScalarBarActor::New();
		scalarBar->SetNumberOfLabels(5);
		// horizontal bar
		//scalarBar->SetOrientationToHorizontal();
		//scalarBar->SetWidth(0.8);
		//scalarBar->SetHeight(0.1);
		// vertical bar
		scalarBar->SetOrientationToVertical();
		scalarBar->SetWidth(0.1);
		scalarBar->SetHeight(1.0);
		
		scalarBar->GetLabelTextProperty()->SetColor(0,0,1);
		scalarBar->GetTitleTextProperty()->SetColor(0,0,1);
		scalarBar->SetLookupTable(mapper->GetLookupTable());
		
		vtkCoordinate* barCoord = vtkCoordinate::New(); // colorbar position in window and dimensions
		barCoord = scalarBar->GetPositionCoordinate();
		barCoord->SetCoordinateSystemToNormalizedViewport();
		barCoord->SetValue(0.1,0.05);
		
		scalarBar->SetLookupTable(lut);
		
		_renderer->AddActor(actor);
		_renderer->AddActor2D(scalarBar);
		_renderer->ResetCamera();
		_renderer->GetRenderWindow()->Render();

		// cleaning up
		/*
		polyData->Delete();
		colors->Delete();
		maskPoints->Delete();
		sphere->Delete();
		glyph->Delete();
		lut->Delete();
		mapper->Delete();
		actor->Delete();
		scalarBar->Delete();
		barCoord->Delete();
		*/
	
	}

	// overwritten virtual function
	void RenderingPointWithColorbar::UpdateYoungModulus(){

		vtkPolyData* polyData = vtkPolyData::New();
		polyData->SetPoints(_points);

		double range[2]; // range for visualization
		_colorBar->GetRange(range);

		// add colors to vtkPolyData
		vtkIntArray* colors = vtkIntArray::New();
		colors->DeepCopy(_colorBar);
		colors->SetName("col");
		polyData->GetPointData()->AddArray(colors);

		
		vtkMaskPoints* maskPoints = vtkMaskPoints::New();
		maskPoints->SetInput(polyData);
		maskPoints->SetOnRatio(5); // interval of points to visualize (one every...)
		
		vtkSphereSource* sphere = vtkSphereSource::New();
		sphere->SetRadius(0.5); // sphere diameter
		
		vtkGlyph3D* glyph = vtkGlyph3D::New();
		glyph->SetInput(maskPoints->GetOutput());
		glyph->SetSource(sphere->GetOutput());

		// lookup table (converts numbers to colors)
		vtkLookupTable* lut = vtkLookupTable::New();
		lut->SetTableRange (range[0], range[1]);
		lut->SetHueRange (0.67, 0.0);
		//lut->SetHueRange (0.0, 0.67);
		lut->Build();
		lut->SetNumberOfColors(256);

		vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
		mapper->SetInput(glyph->GetOutput());
		mapper->SetLookupTable(lut);
		mapper->SetScalarRange(lut->GetTableRange());
		mapper->ScalarVisibilityOn();
		mapper->SetScalarModeToUsePointFieldData(); // SetScalarModeToUsePointData() doesnt allow the points to be colored!!!
		mapper->ColorByArrayComponent("col",0);
		mapper->SetColorModeToMapScalars();
		
		vtkActor* actor = vtkActor::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(_color);

		// colorbar
		vtkScalarBarActor* scalarBar = vtkScalarBarActor::New();
		scalarBar->SetTitle("E [GPa]");
		scalarBar->SetNumberOfLabels(5);
		// horizontal bar
		scalarBar->SetOrientationToHorizontal();
		scalarBar->SetWidth(0.8);
		scalarBar->SetHeight(0.1);
		// vertical bar
		//scalarBar->SetOrientationToVertical();
		//scalarBar->SetWidth(0.1);
		//scalarBar->SetHeight(1.0);
		
		scalarBar->GetLabelTextProperty()->SetColor(0,0,1);
		scalarBar->GetTitleTextProperty()->SetColor(0,0,1);
		scalarBar->SetLookupTable(mapper->GetLookupTable());
		
		vtkCoordinate* barCoord = vtkCoordinate::New(); // colorbar position in window and dimensions
		barCoord = scalarBar->GetPositionCoordinate();
		barCoord->SetCoordinateSystemToNormalizedViewport();
		barCoord->SetValue(0.1,0.05);
		
		scalarBar->SetLookupTable(lut);
		
		_renderer->AddActor(actor);
		_renderer->AddActor2D(scalarBar);
		_renderer->ResetCamera();
		_renderer->GetRenderWindow()->Render();

		// cleaning up
		/*
		polyData->Delete();
		colors->Delete();
		maskPoints->Delete();
		sphere->Delete();
		glyph->Delete();
		lut->Delete();
		mapper->Delete();
		actor->Delete();
		scalarBar->Delete();
		barCoord->Delete();
		*/
	
	}
}
