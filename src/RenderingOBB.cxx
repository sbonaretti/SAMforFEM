/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RenderingOBB.h>

#include <PointWriterXyz.h>

#include <vtkActor.h>
#include <vtkLineSource.h>
#include <vtkOBBTree.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h> //
#include <vtkProperty.h>
#include <vtkRenderWindow.h>

using namespace points;

namespace rendering{
	
	// constructor
	RenderingOBB::RenderingOBB(){

		_mesh = vtkPolyData::New();
		_points = vtkPoints::New();
		
		_color[0] = 1.0; //default color
		_color[1] = 0.0;
		_color[2] = 0.0;

		_one[0] = 0.0; _one[1] = 0.0; _one[2] = 0.0;
		_two[0] = 0.0; _two[1] = 0.0; _two[2] = 0.0;
		_three[0] = 0.0; _three[1] = 0.0; _three[2] = 0.0;
		_four[0] = 0.0; _four[1] = 0.0; _four[2] = 0.0;
		_five[0] = 0.0; _five[1] = 0.0; _five[2] = 0.0;
		_six[0] = 0.0; _six[1] = 0.0; _six[2] = 0.0;
		_seven[0] = 0.0; _seven[1] = 0.0; _seven[2] = 0.0;
		_eight[0] = 0.0; _eight[1] = 0.0; _eight[2] = 0.0;
	}

	
	// destructor
	RenderingOBB::~RenderingOBB(){
	}
	
	// accessors
	void RenderingOBB::GetOBBverteces(double one[], double two[], double three[], double four[], 
									  double five[], double six[], double seven[], double eight[] ){
		
		one[0] = _one[0]; one[1] = _one[1]; one[2] = _one[2];
		two[0] = _two[0]; two[1] = _two[1]; two[2] = _two[2];
		three[0] = _three[0]; three[1] = _three[1]; three[2] = _three[2];
		four[0] = _four[0]; four[1] = _four[1]; four[2] = _four[2];
		five[0] = _five[0]; five[1] = _five[1]; five[2] = _five[2];
		six[0] = _six[0]; six[1] = _six[1]; six[2] = _six[2];
		seven[0] = _seven[0]; seven[1] = _seven[1]; seven[2] = _seven[2];
		eight[0] = _eight[0]; eight[1] = _eight[1]; eight[2] = _eight[2];
	}

	void RenderingOBB::GetOBBcenter(double obbCenter[]){
		
		obbCenter[0] = _obbCenter[0];
		obbCenter[1] = _obbCenter[1];
		obbCenter[2] = _obbCenter[2];
	}

	void RenderingOBB::GetMeshCenterOfMass(double meshCenter[]){
		
		meshCenter[0] = _meshCenterOfMass[0];
		meshCenter[1] = _meshCenterOfMass[1];
		meshCenter[2] = _meshCenterOfMass[2];

	}
	
	// member functions
	void RenderingOBB::CalculateOBB(){

		vtkOBBTree* obb = vtkOBBTree::New();
		double corner[3]; double max[3]; double mid[3]; double min[3]; double size[3];
		if (_mesh->GetNumberOfPoints() != 0)
			_points = _mesh->GetPoints();
		
		obb->ComputeOBB(_points, corner, max, mid, min, size); 

		// OBB vertices
		// corner = _one; max = _two; mid = _three; min = _four;
		_one[0] = corner[0]; _one[1] = corner[1]; _one[2] = corner[2];
		_two[0] = corner[0]+max[0]; _two[1] = corner[1]+max[1]; _two[2] = corner[2]+max[2];
		_three[0] = corner[0]+mid[0]; _three[1] = corner[1]+mid[1]; _three[2] = corner[2]+mid[2];
		_four[0] = corner[0]+min[0]; _four[1] = corner[1]+min[1]; _four[2] = corner[2]+min[2];
		_five[0] = corner[0]+min[0]+max[0]; _five[1] = corner[1]+min[1]+max[1]; _five[2] = corner[2]+min[2]+max[2];
		_six[0] = corner[0]+mid[0]+max[0]; _six[1] = corner[1]+mid[1]+max[1]; _six[2] = corner[2]+mid[2]+max[2];
		_seven[0] = corner[0]+mid[0]+min[0]; _seven[1] = corner[1]+mid[1]+min[1]; _seven[2] = corner[2]+mid[2]+min[2];
		_eight[0] = _two[0] + (_six[0]-_two[0]) + (_five[0]-_two[0]); _eight[1] = _two[1] + (_six[1]-_two[1]) + (_five[1]-_two[1]); _eight[2] = _two[2] + (_six[2]-_two[2]) + (_five[2]-_two[2]);
		

		// checking if the OBB follows the right hand rule
		// OBB sides
		
		double a[3]; 
		a[0] = _four[0] - _one[0]; a[1] = _four[1] - _one[1]; a[2] = _four[2] - _one[2];
		double b[3];
		b[0] = _three[0] - _one[0]; b[1] = _three[1] - _one[1]; b[2] = _three[2] - _one[2];
		// normal
		double normal[3];
		normal[0] = (a[1]*b[2] - a[2]*b[1]) + _one[0];
		normal[1] = (a[0]*b[2] - a[2]*b[0]) + _one[1];
		normal[2] = (a[0]*b[1] - a[1]*b[0]) + _one[2];
		
		
		//double normal[3];
		//normal[0] = (_one[1]*_four[2] - _four[2]_one[1]);
		//normal[1] = (_one[0]-b[2] - _one[2]-b[0]);
		//normal[2] = (_one[0]-b[1] - _one[1]-b[0]);

		std::cout << "normal: " << normal[0] << ' ' << normal[1] << ' ' << normal[2] << std::endl; 

		// CONTROLLARE SE STO PRODOTTO E' GIUSTO : SOMMATO A 1 DEVE DARE 4


		// scalar product with c
		//double direction = normal[0]*c[0] + normal[1]*c[1] + normal[2]*c[2];
		//std::cout << "direction: " << direction << std::endl;
	



		// OBB center
		_obbCenter[0]=(min[0]+mid[0]+max[0])/2 + corner[0];
		_obbCenter[1]=(min[1]+mid[1]+max[1])/2 + corner[1];
		_obbCenter[2]=(min[2]+mid[2]+max[2])/2 + corner[2];

		// mesh center of mass
		int nOfPoints = _mesh->GetNumberOfPoints();
		double x = 0.0; double y = 0.0; double z = 0.0;
		for(int i = 0 ; i < nOfPoints; i++){
			double pt[3];
			_mesh->GetPoint(i,pt);	
			x += pt[0]; y+=pt[1]; z += pt[2];
		}
		
		_meshCenterOfMass[0] = x/nOfPoints;
		_meshCenterOfMass[1] = y/nOfPoints;
		_meshCenterOfMass[2] = z/nOfPoints;

		
	}

	void RenderingOBB::CalculateOBBtop(){

		
		// femur head
		double aSideOne[3]; // 3/4 of the lenght that connects 1 and 2
		for (int i=0; i<3; i++)
			aSideOne[i] = (_two[i]-_one[i]) *3.0/4.0 + _one[i];

		double bSideOne[3]; // 3/4 of the lenght that connects 4 and 5
		for (int i=0; i<3; i++)
			bSideOne[i] = (_five[i]-_four[i]) *3.0/4.0 + _four[i];

		double cSideOne[3]; // 3/4 of the lenght that connects 7 and 8
		for (int i=0; i<3; i++)
			cSideOne[i] = (_eight[i]-_seven[i]) *3.0/4.0 + _seven[i];

		double dSideOne[3]; // 3/4 of the lenght that connects 7 and 8
		for (int i=0; i<3; i++)
			dSideOne[i] = (_six[i]-_three[i]) *3.0/4.0 + _three[i];

		double min[3]; 
		min[0] = aSideOne[0]; 
		if (bSideOne[0] < min[0]) min[0] = bSideOne[0];
		if (cSideOne[0] < min[0]) min[0] = cSideOne[0];
		if (dSideOne[0] < min[0]) min[0] = dSideOne[0];

		min[1] = aSideOne[1]; 
		if (bSideOne[1] < min[1]) min[1] = bSideOne[1]; 
		if (cSideOne[1] < min[1]) min[1] = cSideOne[1]; 
		if (dSideOne[1] < min[1]) min[1] = dSideOne[1];

		min[2] = aSideOne[2]; 
		if (bSideOne[2] < min[2]) min[2] = bSideOne[2];
		if (cSideOne[2] < min[2]) min[2] = cSideOne[2];
		if (dSideOne[2] < min[2]) min[2] = dSideOne[2];

		/*
		double max[3]; 
		max[0] = _two[0];
		if (_five[0] > max[0]) max[0] = _five[0];
		if (_six[0] > max[0]) max[0] = _six[0];
		if (_eight[0] > max[0]) max[0] = _eight[0];

		max[1] = _two[1];
		if (_five[1] > max[1]) max[1] = _five[1];
		if (_six[1] > max[1]) max[1] = _six[1];
		if (_eight[1] > max[1]) max[1] = _eight[1];
		
		max[2] = _two[2];
		if (_five[2] > max[2]) max[2] = _five[2];
		if (_six[2] > max[2]) max[2] = _six[2];
		if (_eight[2] > max[2]) max[2] = _eight[2];
		*/


		// cut the mesh
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			
			double pt[3];
			_mesh->GetPoint(i, pt);

			if (pt[0] > min[0] && pt[1] > min[1] && pt[2] > min[2])
				_points->InsertNextPoint(pt);
		}

		PointWriterXyz* writer = new PointWriterXyz;
		writer->SetInput(_points);
		writer->SetFileName("points.txt");
		writer->Update();

		
		/*
		// vtk obb tree
		// divides the original obb in 2/4/8/... parts in the vertical direction
		// then calculates the nodes that are in that part of the mesh and 
		// finally calculates the obb of those points, without taking into account the original obb
		
		vtkOBBTree* obb = vtkOBBTree::New();
		obb->SetDataSet(_mesh);
		obb->BuildLocator();
		vtkPolyData* obbBox = vtkPolyData::New();
		obb->GenerateRepresentation(3,obbBox);
		vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
		writer->SetInput(obbBox);
		writer->SetFileName("obbBoxTest.vtk");
		writer->Update();
		*/

	}

	void RenderingOBB::CalculateOBBbottom(){

		
		// femur head
		double aSideOne[3]; // 1/4 of the lenght that connects 1 and 2
		for (int i=0; i<3; i++)
			aSideOne[i] = (_two[i]-_one[i]) *1.0/4.0 + _one[i];

		double bSideOne[3]; // 1/4 of the lenght that connects 4 and 5
		for (int i=0; i<3; i++)
			bSideOne[i] = (_five[i]-_four[i]) *1.0/4.0 + _four[i];

		double cSideOne[3]; // 1/4 of the lenght that connects 7 and 8
		for (int i=0; i<3; i++)
			cSideOne[i] = (_eight[i]-_seven[i]) *1.0/4.0 + _seven[i];

		double dSideOne[3]; // 1/4 of the lenght that connects 7 and 8
		for (int i=0; i<3; i++)
			dSideOne[i] = (_six[i]-_three[i]) *1.0/4.0 + _three[i];

		double min[3]; 
		min[0] = aSideOne[0]; 
		if (bSideOne[0] < min[0]) min[0] = bSideOne[0];
		if (cSideOne[0] < min[0]) min[0] = cSideOne[0];
		if (dSideOne[0] < min[0]) min[0] = dSideOne[0];

		min[1] = aSideOne[1]; 
		if (bSideOne[1] < min[1]) min[1] = bSideOne[1]; 
		if (cSideOne[1] < min[1]) min[1] = cSideOne[1]; 
		if (dSideOne[1] < min[1]) min[1] = dSideOne[1];

		min[2] = aSideOne[2]; 
		if (bSideOne[2] < min[2]) min[2] = bSideOne[2];
		if (cSideOne[2] < min[2]) min[2] = cSideOne[2];
		if (dSideOne[2] < min[2]) min[2] = dSideOne[2];
		

		std::cout << aSideOne[0] << ' ' << aSideOne[1] << ' ' << aSideOne[2] << std::endl;
		std::cout << bSideOne[0] << ' ' << bSideOne[1] << ' ' << bSideOne[2] << std::endl;
		std::cout << cSideOne[0] << ' ' << cSideOne[1] << ' ' << cSideOne[2] << std::endl;
		std::cout << dSideOne[0] << ' ' << dSideOne[1] << ' ' << dSideOne[2] << std::endl;


		std::cout << min[0] << ' ' << min[1] << ' ' << min[2] << std::endl;




		/*
		double max[3]; 
		max[0] = _two[0];
		if (_five[0] > max[0]) max[0] = _five[0];
		if (_six[0] > max[0]) max[0] = _six[0];
		if (_eight[0] > max[0]) max[0] = _eight[0];

		max[1] = _two[1];
		if (_five[1] > max[1]) max[1] = _five[1];
		if (_six[1] > max[1]) max[1] = _six[1];
		if (_eight[1] > max[1]) max[1] = _eight[1];
		
		max[2] = _two[2];
		if (_five[2] > max[2]) max[2] = _five[2];
		if (_six[2] > max[2]) max[2] = _six[2];
		if (_eight[2] > max[2]) max[2] = _eight[2];
		*/


		// cut the mesh
		_points = vtkPoints::New();
		for (int i=0; i<_mesh->GetNumberOfPoints(); i++){
			
			double pt[3];
			_mesh->GetPoint(i, pt);

			if (pt[0] < min[0] && pt[1] < min[1] && pt[2] < min[2])
				_points->InsertNextPoint(pt);
		}

		PointWriterXyz* writer = new PointWriterXyz;
		writer->SetInput(_points);
		writer->SetFileName("points.txt");
		writer->Update();
	}

	
	// overwritten virtual function (obb rendering)
	void RenderingOBB::Update(){

		
		double lineWidth = 0.5;
		
		// one and two
		vtkLineSource* line1 = vtkLineSource::New();
		line1->SetPoint1(_one);
		line1->SetPoint2(_two);
		vtkPolyDataMapper * mapper1 = vtkPolyDataMapper ::New();
		mapper1->SetInput(line1->GetOutput());
		vtkActor * actor1 = vtkActor::New();
		actor1->SetMapper(mapper1);
		_color[0] = 0.0; _color[1] = 0.0; _color[2] = 1.0;
		actor1->GetProperty()->SetColor(_color);
		actor1->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor1);

		// one and three
		vtkLineSource* line2 = vtkLineSource::New();
		line2->SetPoint1(_one);
		line2->SetPoint2(_three);
		vtkPolyDataMapper * mapper2 = vtkPolyDataMapper ::New();
		mapper2->SetInput(line2->GetOutput());
		vtkActor * actor2 = vtkActor::New();
		actor2->SetMapper(mapper2);
		_color[0] = 0.0; _color[1] = 1.0; _color[2] = 0.0;
		actor2->GetProperty()->SetColor(_color);
		actor2->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor2);

		// one and four
		vtkLineSource* line3 = vtkLineSource::New();
		line3->SetPoint1(_one);
		line3->SetPoint2(_four);
		vtkPolyDataMapper * mapper3 = vtkPolyDataMapper ::New();
		mapper3->SetInput(line3->GetOutput());
		vtkActor * actor3 = vtkActor::New();
		actor3->SetMapper(mapper3);
		_color[0] = 1.0; _color[1] = 0.0; _color[2] = 0.0;
		actor3->GetProperty()->SetColor(_color);
		actor3->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor3);

		// four and five
		vtkLineSource* line4 = vtkLineSource::New();
		line4->SetPoint1(_four);
		line4->SetPoint2(_five);
		vtkPolyDataMapper * mapper4 = vtkPolyDataMapper ::New();
		mapper4->SetInput(line4->GetOutput());
		vtkActor * actor4 = vtkActor::New();
		actor4->SetMapper(mapper4);
		_color[0] = 0.0; _color[1] = 0.0; _color[2] = 0.0;
		actor4->GetProperty()->SetColor(_color);
		actor4->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor4);

		// two and five
		vtkLineSource* line5 = vtkLineSource::New();
		line5->SetPoint1(_two);
		line5->SetPoint2(_five);
		vtkPolyDataMapper * mapper5 = vtkPolyDataMapper ::New();
		mapper5->SetInput(line5->GetOutput());
		vtkActor * actor5 = vtkActor::New();
		actor5->SetMapper(mapper5);
		actor5->GetProperty()->SetColor(_color);
		actor5->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor5);

		// two and six
		vtkLineSource* line6 = vtkLineSource::New();
		line6->SetPoint1(_two);
		line6->SetPoint2(_six);
		vtkPolyDataMapper * mapper6 = vtkPolyDataMapper ::New();
		mapper6->SetInput(line6->GetOutput());
		vtkActor * actor6 = vtkActor::New();
		actor6->SetMapper(mapper6);
		actor6->GetProperty()->SetColor(_color);
		actor6->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor6);

		// three and six
		vtkLineSource* line7 = vtkLineSource::New();
		line7->SetPoint1(_three);
		line7->SetPoint2(_six);
		vtkPolyDataMapper * mapper7 = vtkPolyDataMapper ::New();
		mapper7->SetInput(line7->GetOutput());
		vtkActor * actor7 = vtkActor::New();
		actor7->SetMapper(mapper7);
		actor7->GetProperty()->SetColor(_color);
		actor7->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor7);

		// three and seven
		vtkLineSource* line8 = vtkLineSource::New();
		line8->SetPoint1(_three);
		line8->SetPoint2(_seven);
		vtkPolyDataMapper * mapper8 = vtkPolyDataMapper ::New();
		mapper8->SetInput(line8->GetOutput());
		vtkActor * actor8 = vtkActor::New();
		actor8->SetMapper(mapper8);
		actor8->GetProperty()->SetColor(_color);
		actor8->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor8);

		// four and seven
		vtkLineSource* line9 = vtkLineSource::New();
		line9->SetPoint1(_four);
		line9->SetPoint2(_seven);
		vtkPolyDataMapper * mapper9 = vtkPolyDataMapper ::New();
		mapper9->SetInput(line9->GetOutput());
		vtkActor * actor9 = vtkActor::New();
		actor9->SetMapper(mapper9);
		actor9->GetProperty()->SetColor(_color);
		actor9->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor9);

		// seven and eight
		vtkLineSource* line10 = vtkLineSource::New();
		line10->SetPoint1(_seven);
		line10->SetPoint2(_eight);
		vtkPolyDataMapper * mapper10 = vtkPolyDataMapper ::New();
		mapper10->SetInput(line10->GetOutput());
		vtkActor * actor10 = vtkActor::New();
		actor10->SetMapper(mapper10);
		actor10->GetProperty()->SetColor(_color);
		actor10->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor10);

		// six and eight
		vtkLineSource* line11 = vtkLineSource::New();
		line11->SetPoint1(_six);
		line11->SetPoint2(_eight);
		vtkPolyDataMapper * mapper11 = vtkPolyDataMapper ::New();
		mapper11->SetInput(line11->GetOutput());
		vtkActor * actor11 = vtkActor::New();
		actor11->SetMapper(mapper11);
		actor11->GetProperty()->SetColor(_color);
		actor11->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor11);

		// six and eight
		vtkLineSource* line12 = vtkLineSource::New();
		line12->SetPoint1(_five);
		line12->SetPoint2(_eight);
		vtkPolyDataMapper * mapper12 = vtkPolyDataMapper ::New();
		mapper12->SetInput(line12->GetOutput());
		vtkActor * actor12 = vtkActor::New();
		actor12->SetMapper(mapper12);
		actor12->GetProperty()->SetColor(_color);
		actor12->GetProperty()->SetLineWidth(lineWidth);
		_renderer->AddActor(actor12);
	}
	
}