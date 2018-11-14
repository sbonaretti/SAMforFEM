/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshSimplifyAndSmooth.h>

#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>

#include <QString>
#include <QFile>
#include <QTextStream>
#include <QProcess>

namespace mesh{

	// constructor
	MeshSimplifyAndSmooth::MeshSimplifyAndSmooth(){

		_simplifyValue = 0.0;
		_smoothNofIteration = 0;
	}

	// destructor
	MeshSimplifyAndSmooth::~MeshSimplifyAndSmooth(){
	}

	// member functions
	void MeshSimplifyAndSmooth::Update(){

		std::cout << "simplifyValue: " << _simplifyValue << std::endl;
		std::cout << "smoothNofIteration: " << _smoothNofIteration<< std::endl;

		// fileNames and path
		QString csvFileName = ("C:/Program Files (x86)/MRFSurface/nodes.csv");
		QString MRFSurfacePath = ("C:/Program Files (x86)/MRFSurface/MRFSurface.exe");
		QString meshFileName = ("C:/Program Files (x86)/MRFSurface/mesh.stl");
		int waitForMRFSurfaceFinished = 2000000;
			
		// original mesh
		std::cout << "number of nodes of the original mesh: " <<  _mesh->GetNumberOfPoints() << std::endl;
		vtkSTLWriter* meshWriter = vtkSTLWriter::New();
		meshWriter->SetFileName("original.stl");
		meshWriter->SetInput(_mesh);
		//meshWriter->Update();
		
		
		// simplify mesh
		vtkDecimatePro* quadricDecimation = vtkDecimatePro::New();
		if (_simplifyValue != 0.0){
			std::cout << "simplify mesh" << std::endl;	
			
			quadricDecimation->SetInput(_mesh);
			quadricDecimation->SetTargetReduction(_simplifyValue);
			quadricDecimation->Update();
			_mesh = quadricDecimation->GetOutput();
			meshWriter->SetFileName("simplified.stl");
			meshWriter->SetInput(quadricDecimation->GetOutput());
			//meshWriter->Update();
			std::cout << "number of nodes of the decimated mesh: " <<  quadricDecimation->GetOutput()->GetNumberOfPoints() << std::endl;
		}
			
		
		// smooth mesh
		vtkSmoothPolyDataFilter* smoother = vtkSmoothPolyDataFilter::New(); 
		if (_smoothNofIteration != 0){
			std::cout << "smooth mesh" << std::endl;	
			smoother->SetInput(_mesh);					
			smoother->SetNumberOfIterations(_smoothNofIteration);
			smoother->Update();
			_mesh = smoother->GetOutput();
			meshWriter->SetFileName("smoothed.stl");
			meshWriter->SetInput(_mesh);
			//meshWriter->Update();
		}

	
		// recreate topology
		if ((_simplifyValue != 0.0) || (_smoothNofIteration != 0)){
	
			// write nodes coordinates (input for MRFSurface)
			QFile outFile("test.csv");
			//QFile outFile(csvFileName.toAscii().data());
			outFile.open(QIODevice::WriteOnly | QIODevice::Text);
			QTextStream file(&outFile);
			//std::cout << _mesh->GetNumberOfPoints() << std::endl;
			for (int a=0; a<_mesh->GetNumberOfPoints(); a++){
				double pt[3];
				_mesh->GetPoint(a, pt);
				file << pt[0] << ',' << pt[1] << ',' << pt[2] << endl;
			}
			outFile.close();

			// execute MRFSurface
			std::cout << "run MRFSurface" << std::endl;
			
			QStringList arguments;
			arguments << "-t" << "0" << "-i" << "test.csv" << "-o" << "mesh.stl"; // 0 corresponds to the point cloud
			//arguments << "-t" << "0" << "-i" << csvFileName.toAscii().data() << "-o" << meshFileName.toAscii().data(); // 0 corresponds to the point cloud
			
			QProcess myProcess;
			myProcess.start(MRFSurfacePath, arguments);
			myProcess.waitForFinished(waitForMRFSurfaceFinished);
			
			// cell inversion
			vtkSTLReader* meshReader = vtkSTLReader::New();
			meshReader->SetFileName("mesh.stl");
			//meshReader->SetFileName(meshFileName.toAscii().data());
			meshReader->Update();
			for (int a=0; a<meshReader->GetOutput()->GetNumberOfCells(); a++)
				meshReader->GetOutput()->ReverseCell(a);
			_mesh = meshReader->GetOutput();

		}

		//cleaning up
		quadricDecimation->Delete();
		smoother->Delete();
		//meshWriter->Delete(); //it crashes
		//meshReader->Delete(); //it crashes

	
	
	}

}