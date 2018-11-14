/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <RegistrationValidationSurfaceMesh.h>

#include <MeshExtractOuterSurface.h>
#include <MeshReaderMorpherVolume.h>
#include <PointWriterXyz.h>

#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>

#include <QFile>
#include <QTextStream>

using namespace mesh;
using namespace points;


namespace registration{

	// constructor
	RegistrationValidationSurfaceMesh::RegistrationValidationSurfaceMesh(){

		_flag = 0;
	}

	// destructor
	RegistrationValidationSurfaceMesh::~RegistrationValidationSurfaceMesh(){
	}

	// member function
	void RegistrationValidationSurfaceMesh::Update(){

		vtkPolyData* registeredMesh = vtkPolyData::New();
		vtkSTLReader* stlReader = vtkSTLReader::New();

		for (int i=0; i<_registeredFileNames.size(); i++){

			
			if (_flag == 0){ // for the mesh tool
				
				// load the registered mesh (volume .cdb)
				MeshReaderMorpherVolume* meshReaderMorpherVolume = new MeshReaderMorpherVolume;
				meshReaderMorpherVolume->SetFileName(_registeredFileNames[i].ascii());
				std::cout << std::endl;
				std::cout << "registered mesh " << i+1 << ": " << _registeredFileNames[i].ascii() << std::endl;
				meshReaderMorpherVolume->Update();

				// extract the outer surface
				std::cout << "extract outer nodes from the registered mesh" << std::endl;
				MeshExtractOuterSurface* extractOuterSurface = new MeshExtractOuterSurface;
				extractOuterSurface->SetVolumeMesh(meshReaderMorpherVolume->GetOutput());
				extractOuterSurface->Update();
				registeredMesh->DeepCopy(extractOuterSurface->GetSurfaceMesh());

				delete meshReaderMorpherVolume;
				delete extractOuterSurface;
			}

			else if (_flag == 1){ // for the image tool
				
				//load the registered mesh (surface .stl)
				std::cout << std::endl;
				std::cout << "registered mesh " << i+1 << ": " << _registeredFileNames[i].ascii() << std::endl;
				stlReader->SetFileName(_registeredFileNames[i].ascii());
				stlReader->Update();
				std::cout << "number of nodes: " << stlReader->GetOutput()->GetNumberOfPoints() << std::endl;
				std::cout << "number of elements: " << stlReader->GetOutput()->GetNumberOfCells() << std::endl;
				registeredMesh->DeepCopy(stlReader->GetOutput());
			}

			
			// load the correspondent original .stl
			QString temp = _registeredFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("\\")+1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("/")+1);
			}
			for (int y=0; y<_originalFileNames.size(); y++){
				if (_originalFileNames[y].contains(temp)){
					std::cout << "original mesh: " << _originalFileNames[y].ascii() << std::endl;
					// load the mesh
					stlReader->SetFileName(_originalFileNames[y].ascii());
					stlReader->Update();
					std::cout << "number of nodes: " << stlReader->GetOutput()->GetNumberOfPoints() << std::endl;
					std::cout << "number of elements: " << stlReader->GetOutput()->GetNumberOfCells() << std::endl;
				}
			}

			// find closest points for the point-to-surface difference: extracted mesh -> original stl
			std::cout << "number of points: " << registeredMesh->GetNumberOfPoints() << std::endl;
			double extractedPointsToOriginalStlMaxDistance;
			findPointToSurfaceDistances(registeredMesh->GetPoints(), stlReader->GetOutput(), extractedPointsToOriginalStlMaxDistance);
			std::cout << "maximum distance registered mesh -> original mesh: " << extractedPointsToOriginalStlMaxDistance << std::endl;

			// find closest points for the point-to-surface difference: original stl -> extracted mesh 
			std::cout << "number of points: " << stlReader->GetOutput()->GetPoints()->GetNumberOfPoints() << std::endl;
			double originalPointsToExtractedStlMaxDistance;
			findPointToSurfaceDistances(stlReader->GetOutput()->GetPoints(), registeredMesh, originalPointsToExtractedStlMaxDistance);
			std::cout << "maximum distance original mesh -> registered mesh: " << originalPointsToExtractedStlMaxDistance << std::endl;

			// set the Hausdorff distance in the array
			if (extractedPointsToOriginalStlMaxDistance > originalPointsToExtractedStlMaxDistance)
				_HausdorffDistances->InsertNextValue(extractedPointsToOriginalStlMaxDistance);
			else
				_HausdorffDistances->InsertNextValue(originalPointsToExtractedStlMaxDistance);

			// cleaning up
			//stlReader->Delete();
			

		}

		// calculate the average of the Hausdorff distances
		for (int i=0; i<_HausdorffDistances->GetNumberOfTuples(); i++) {
			_average += _HausdorffDistances->GetValue(i);
		}
		_average /= _HausdorffDistances->GetNumberOfTuples();
		std::cout << std::endl;
		std::cout << "average: " << _average << std::endl;
		
		// standard deviation
		for (int i=0; i<_HausdorffDistances->GetNumberOfTuples(); i++)
			_standardDeviation += (_HausdorffDistances->GetValue(i) - _average)*(_HausdorffDistances->GetValue(i) - _average);
		_standardDeviation /= _HausdorffDistances->GetNumberOfTuples()-1;
		_standardDeviation = sqrt (_standardDeviation);
		std::cout << "standardDeviation: " << _standardDeviation << std::endl;

		// standard error
		_standardError = _standardDeviation / sqrt(double(_HausdorffDistances->GetNumberOfTuples()));
		std::cout << "standardError: " << _standardError << std::endl;

		// write it to file
		QFile outFile("hausdorff distance mesh.txt");
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&outFile);
		writeFile << "average: " << _average << endl;
		writeFile << "standardDeviation: " << _standardDeviation << endl;
		writeFile << "standardError: " << _standardError << endl;
		outFile.close();


	
	}
}