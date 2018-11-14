/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RegistrationProcrustesAlignment.h>

#include <MeshReaderMorpherVolume.h>
#include <MeshWriterAnsys.h>

#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkProcrustesAlignmentFilter.h>
#include <vtkTransformPolyDataFilter.h>


using namespace mesh;

namespace registration{

	// constructor
	RegistrationProcrustesAlignment::RegistrationProcrustesAlignment(){
		
		_flag = 0;
	
	}

	// destructor
	RegistrationProcrustesAlignment::~RegistrationProcrustesAlignment(){
	}

	// member function
	void RegistrationProcrustesAlignment::Update(){
		
		vtkPolyData* topologySaver = vtkPolyData::New();
		
		vtkProcrustesAlignmentFilter* procrustes = vtkProcrustesAlignmentFilter::New();
		procrustes->SetNumberOfInputs(_fileNames.size());
		if (_flag == 0){
			std::cout << "procrustes - rigid alignment" << std::endl;
			procrustes->GetLandmarkTransform()->SetModeToRigidBody();
		}
		else if (_flag == 1){
			std::cout << "procrustes - affine alignment" << std::endl;
			procrustes->GetLandmarkTransform()->SetModeToAffine();
		}
		else 
			std::cout << "procrustes - similarity alignment" << std::endl;
		
		std::cout << "mesh to align:" << std::endl;
		for (int i=0; i<_fileNames.size(); i++){
			
			// read meshes
			std::cout << _fileNames[i].toAscii().data() << std::endl;
			MeshReaderMorpherVolume* reader = new MeshReaderMorpherVolume;
			reader->SetFileName(_fileNames[i]);
			reader->Update();
			
			/*
			// find the mesh center of mass
			double centerOfMass[3]; centerOfMass[0]=0.0; centerOfMass[1]=0.0; centerOfMass[2]=0.0;
			for (int y=0; y<reader->GetOutput()->GetNumberOfPoints(); y++){
				double point[3];
				reader->GetOutput()->GetPoint(y, point);
				centerOfMass[0] += point[0];
				centerOfMass[1] += point[1];
				centerOfMass[2] += point[2];
			}
			centerOfMass[0] /= reader->GetOutput()->GetNumberOfPoints();
			centerOfMass[1] /= reader->GetOutput()->GetNumberOfPoints();
			centerOfMass[2] /= reader->GetOutput()->GetNumberOfPoints();
		
			// move the mesh center of mass to the origin (0,0,0)
			vtkPoints* points = vtkPoints::New();		
			for (int y=0; y<reader->GetOutput()->GetNumberOfPoints(); y++){
				double point[3];
				reader->GetOutput()->GetPoint(y, point);
				point[0] -= centerOfMass[0];
				point[1] -= centerOfMass[1];
				point[2] -= centerOfMass[2];
				points->InsertNextPoint(point);
			}
			reader->GetOutput()->SetPoints(points);
			*/
						
			// give it to filter
			procrustes->SetInput(i, reader->GetOutput());
			// saving topology
			if (i == 0)
				topologySaver = reader->GetOutput();
			
			delete reader;
		}
		
		// align
		std::cout << "run procrustes alignment" << std::endl;
		procrustes->Update();
			
		// save aligned meshes
		std::cout << "save meshes" << std::endl;
		for (int i=0; i<_fileNames.size(); i++){
			
			_fileNames[i].replace(QString (".cdb"), QString ("_aligned.cdb"));
			std::cout << _fileNames[i].toAscii().data() << std::endl;
			
			MeshWriterAnsys* writer = new MeshWriterAnsys;
			writer->SetFileName(_fileNames[i]);
			writer->MeshOn();
			topologySaver->SetPoints(procrustes->GetOutput(i)->GetPoints());
			writer->SetMesh(topologySaver);
			writer->Update();
			
			delete writer;
		}

	
	}

	void RegistrationProcrustesAlignment::IterativeClosestPoints(){
	
		vtkIterativeClosestPointTransform* icp = vtkIterativeClosestPointTransform::New();
		icp->SetSource(_sourceMesh);
		icp->SetTarget(_targetMesh);
		icp->GetLandmarkTransform()->SetModeToRigidBody();
		icp->SetMaximumNumberOfIterations(50);
		icp->StartByMatchingCentroidsOn();
		icp->Modified();
		icp->Update();
 
		// Get the resulting transformation matrix (this matrix takes the source points to the target points)
		vtkMatrix4x4* m = icp->GetMatrix();
		std::cout << "The resulting matrix is: " << *m << std::endl;
 
		// Transform the source points by the ICP solution
		vtkTransformPolyDataFilter* icpTransformFilter = vtkTransformPolyDataFilter::New();
		icpTransformFilter->SetInput(_sourceMesh);
		icpTransformFilter->SetTransform(icp);
		icpTransformFilter->Update();

		_outputMesh = icpTransformFilter->GetOutput();




	}

}