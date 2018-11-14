/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <ImageHandler.h>
#include <MeshSimplifyAndSmooth.h>
#include <RegistrationValidationSurfaceImage.h>

#include <vtkSTLWriter.h>

#include <QFile>
#include <QTextStream>

using namespace image;
using namespace mesh;


namespace registration{

	// constructor
	RegistrationValidationSurfaceImage::RegistrationValidationSurfaceImage(){
	}

	// destructor
	RegistrationValidationSurfaceImage::~RegistrationValidationSurfaceImage(){
	}

	// member function
	void RegistrationValidationSurfaceImage::Update(){

		vtkSTLWriter* meshWriter = vtkSTLWriter::New();
		
		for (int i=0; i<_registeredFileNames.size(); i++){
			
			
			// load the registered image
			ImageHandler* imageHandler = new ImageHandler;
			std::cout << std::endl;
			std::cout << "registered image: " << _registeredFileNames[i].ascii() << std::endl;
			imageHandler->SetImageFileName(_registeredFileNames[i].ascii());
			imageHandler->MetafileReaderUpdate();
			
			// extract mesh, smooth and simplify it
			imageHandler->MarchingCubesUpdate();
			MeshSimplifyAndSmooth* registeredMesh = new MeshSimplifyAndSmooth;
			registeredMesh->SetInput(imageHandler->GetMesh());
			registeredMesh->SetSimplifyValue(0.4);
			registeredMesh->SetSmoothIteration(150);
			registeredMesh->Update();
			// save mesh
			_registeredFileNames[i].replace(QString(".mhd"), QString(".stl"));
			meshWriter->SetInput(registeredMesh->GetOutput());
			meshWriter->SetFileName(_registeredFileNames[i].toAscii());
			//meshWriter->Update();

			
			// load the correspondent original image
			QString temp = _registeredFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("\\")+1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("/")+1);
			}
			QString meshTemp;
			for (int y=0; y<_originalFileNames.size(); y++){
				if (_originalFileNames[y].contains(temp)){
					std::cout << "original image: " << _originalFileNames[y].ascii() << std::endl;
					// load the mesh
					imageHandler->SetImageFileName(_originalFileNames[y].ascii());
					imageHandler->MetafileReaderUpdate();
					_originalFileNames[y].replace(QString(".mhd"), QString(".stl")); // to save the mesh below
					meshTemp = _originalFileNames[y];
				}
			}

			// extract mesh, smooth and simplify it
			imageHandler->MarchingCubesUpdate();
			MeshSimplifyAndSmooth* originalMesh = new MeshSimplifyAndSmooth;
			originalMesh->SetInput(imageHandler->GetMesh());
			originalMesh->SetSimplifyValue(0.4);
			originalMesh->SetSmoothIteration(150);
			originalMesh->Update();
			// save mesh
			meshWriter->SetInput(originalMesh->GetOutput());
			meshWriter->SetFileName(meshTemp.toAscii().data());
			//meshWriter->Update();

			// find closest points for the point-to-surface difference: registered mesh -> original mesh
			std::cout << "number of points: " << registeredMesh->GetOutput()->GetPoints()->GetNumberOfPoints() << std::endl;
			double extractedPointsToOriginalStlMaxDistance;
			findPointToSurfaceDistances(registeredMesh->GetOutput()->GetPoints(), originalMesh->GetOutput(), extractedPointsToOriginalStlMaxDistance);
			std::cout << "maximum distance registered mesh -> original mesh: " << extractedPointsToOriginalStlMaxDistance << std::endl;

			// find closest points for the point-to-surface difference: original mesh -> registered mesh 
			std::cout << "number of points: " << originalMesh->GetOutput()->GetPoints()->GetNumberOfPoints() << std::endl;
			double originalPointsToExtractedStlMaxDistance;
			findPointToSurfaceDistances(originalMesh->GetOutput()->GetPoints(), registeredMesh->GetOutput(), originalPointsToExtractedStlMaxDistance);
			std::cout << "maximum distance original mesh -> registered mesh: " << originalPointsToExtractedStlMaxDistance << std::endl;

			// set the Hausdorff distance in the array
			if (extractedPointsToOriginalStlMaxDistance > originalPointsToExtractedStlMaxDistance)
				_HausdorffDistances->InsertNextValue(extractedPointsToOriginalStlMaxDistance);
			else
				_HausdorffDistances->InsertNextValue(originalPointsToExtractedStlMaxDistance);
			
			// cleaning up
			delete imageHandler;
			delete registeredMesh;
			delete originalMesh;
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
		QFile outFile("hausdorff distance image.txt");
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFile(&outFile);
		writeFile << "average: " << _average << endl;
		writeFile << "standardDeviation: " << _standardDeviation << endl;
		writeFile << "standardError: " << _standardError << endl;
		outFile.close();

		// cleaning up
		meshWriter->Delete();
	}
}

