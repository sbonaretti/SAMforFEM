/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshQuality.h>

#include <MeshReaderAbaqus.h>
#include <MeshReaderMorpherVolume.h>

#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkMeshQuality.h>
#include <vtkTetra.h>

#include <VnlWriterMatrix.h>
#include <VnlWriterVector.h>
#include <vnl/vnl_vector.h>

using namespace vnl;

namespace mesh{

	// constructor
	MeshQuality::MeshQuality(){	
	}

	// destructor
	MeshQuality::~MeshQuality(){
	}

	// member functions
	void MeshQuality::CalculateJacobian(){

		// define if the meshes are from abaqus (.inp) or ansys (.cdb)
		bool fileType;
		if (_meshFileNames[0].endsWith(".inp"))
			fileType = true;
		else if (_meshFileNames[0].endsWith(".cdb"))
			fileType = false;
		else
			std::cout << "input files not supported" << std::endl;

		// folder path
		QString folderName = _meshFileNames[0];
		
		// jacobian
		_jacobiansQuality.set_size(2, _meshFileNames.size());


		// for each mesh
		for (int i=0; i<_meshFileNames.size(); i++){

			vtkPolyData* mesh = vtkPolyData::New();
			
			// if abaqus 
			if (fileType == true){
				// read mesh
				MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
				meshReader->SetFileName(_meshFileNames[i]);
				std::cout << _meshFileNames[i].toAscii().data() << std::endl; 
				meshReader->Update();
				mesh->DeepCopy(meshReader->GetOutput());
				delete meshReader;
			}
			
			// if ansys
			else if (fileType == false){
				// read mesh
				MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
				meshReader->SetFileName(_meshFileNames[i]);
				std::cout << _meshFileNames[i].toAscii().data() << std::endl; 
				meshReader->Update();
				mesh->DeepCopy(meshReader->GetOutput());
				delete meshReader;
			}
			
		
			// extracting cells
			vtkGenericCell* cell = vtkGenericCell::New();
			vtkIdList* idList;
			vtkTetra* tetra = vtkTetra::New();
			vtkMeshQuality* quality = vtkMeshQuality::New();
		
			int count = 0;
			_jacobians.set_size(mesh->GetNumberOfCells());
			
			for (int y=0; y<mesh->GetNumberOfCells(); y++){
				
				// get mesh element
				mesh->GetCell(y, cell);
				idList = cell->GetPointIds();
				
				// put element vertices in the tetra
				double pt0[3]; mesh->GetPoint(idList->GetId(0),pt0); tetra->GetPoints()->SetPoint(0, pt0);
				double pt1[3]; mesh->GetPoint(idList->GetId(1),pt1); tetra->GetPoints()->SetPoint(1, pt1);
				double pt2[3]; mesh->GetPoint(idList->GetId(2),pt2); tetra->GetPoints()->SetPoint(2, pt2);
				double pt3[3]; mesh->GetPoint(idList->GetId(3),pt3); tetra->GetPoints()->SetPoint(3, pt3);

				// put element connections in the tetra
				tetra->GetPointIds()->SetId(0,0);
				tetra->GetPointIds()->SetId(1,1);
				tetra->GetPointIds()->SetId(2,2);
				tetra->GetPointIds()->SetId(3,3);

				// calculate the jacobian
				_jacobians(y) = quality->TetScaledJacobian(tetra);
				if (quality->TetScaledJacobian(tetra) <= 0){
					count ++;	
				}
			}	

			// jacobian file name
			QString fileName = _meshFileNames[i];
			if (fileType == true)
				fileName.replace(QString(".inp"),QString("_jacobian.txt"));
			else if (fileType == false)
				fileName.replace(QString(".cdb"),QString("_jacobian.txt"));
			// writing		
			VnlWriterVector* vectorWriter = new VnlWriterVector;
			vectorWriter->SetVnlVector(_jacobians);
			vectorWriter->SetFileName(fileName);
			vectorWriter->Update();
	
			// extracting mesh id name (original meshes)
			/*QString temp = _meshFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("\\")+1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("/")+1);
			}*/
			// extracting mesh id name (created instances)
			QString temp = _meshFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(0, temp.lastIndexOf("\\")+10);
				temp.remove(".cdb");
			}
			else {
				temp.remove(0, temp.lastIndexOf("/")+10);
				temp.remove(".cdb");		
			}	
			
			// give it to jacobian quality
			_jacobiansQuality(0,i) = temp.toDouble();
			// number of 
			_jacobiansQuality(1,i) = count;

			quality->Delete();
			mesh->Delete();
			delete vectorWriter;
		}
	
		// info about the jacobians
		// number of wrong meshes
		int count = 0;
		for (int i=0; i<_meshFileNames.size(); i++){
			if (_jacobiansQuality(1,i) != 0)
				count++;
		}
		std::cout << "number of meshes with negative or null jacobian elements: " << count << std::endl;
		// create an array containing the null-negative data
		vnl_vector<double> negativeJacobians;
		negativeJacobians.set_size(count);
		count = 0;
		for (int i=0; i<_meshFileNames.size(); i++){
			if (_jacobiansQuality(1,i) != 0){
				negativeJacobians(count) = _jacobiansQuality(1,i);
				count++;
			}
		}
		std::cout << "negative jacobians: " << negativeJacobians << std::endl;
		// average
		std::cout << "average number of null or negative elements: " << negativeJacobians.mean() << std::endl;
		// std dev
		double stdDev = 0.0;
		for (int i=0; i<negativeJacobians.size(); i++)
			stdDev += (negativeJacobians(i) - negativeJacobians.mean())*(negativeJacobians(i) - negativeJacobians.mean());
		stdDev /= (negativeJacobians.size()-1);
		stdDev = sqrt (stdDev);
		std::cout << "standard deviation: " << stdDev << std::endl;
		
		// write matrices and vector
		if (folderName.lastIndexOf("/") == -1){
			folderName.remove(folderName.lastIndexOf("\\")+1,folderName.size()-1);
			}
		else {
			folderName.remove(folderName.lastIndexOf("/")+1,folderName.size()-1);
		}
		folderName.append("JacobiansQuality.txt");
		VnlWriterMatrix* writer = new VnlWriterMatrix; 
		writer->SetVnlMatrix(_jacobiansQuality.transpose());
		writer->SetFileName(folderName);
		writer->MatrixShapeUpdate();

		/*
		VnlWriterVector* vectorWriter = new VnlWriterVector;
		folderName.replace(QString("JacobiansQuality.txt"), QString("NegativeJacobians.txt"));
		vectorWriter->SetVnlVector(negativeJacobians);
		vectorWriter->SetFileName(folderName);
		vectorWriter->Update();
		*/
		delete writer;

		std::cout << "files written" << std::endl;
	}

	void MeshQuality::CalculateEdgeRatio(){

		// define if the meshes are from abaqus (.inp) or ansys (.cdb)
		bool fileType;
		if (_meshFileNames[0].endsWith(".inp"))
			fileType = true;
		else if (_meshFileNames[0].endsWith(".cdb"))
			fileType = false;
		else
			std::cout << "input files not supported" << std::endl;

		// folder path
		QString folderName = _meshFileNames[0];
		
		// edge ratio quality
		_edgeRatioQuality.set_size(3, _meshFileNames.size()); // mesh id, edge ratio average, edge ratio std dev


		// for each mesh
		for (int i=0; i<_meshFileNames.size(); i++){

			vtkPolyData* mesh = vtkPolyData::New();
			
			// if abaqus 
			if (fileType == true){
				// read mesh
				MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
				meshReader->SetFileName(_meshFileNames[i]);
				std::cout << _meshFileNames[i].toAscii().data() << std::endl; 
				meshReader->Update();
				mesh->DeepCopy(meshReader->GetOutput());

				delete meshReader;
			}
			
			// if ansys
			else if (fileType == false){
				// read mesh
				MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
				meshReader->SetFileName(_meshFileNames[i]);
				std::cout << _meshFileNames[i].toAscii().data() << std::endl; 
				meshReader->Update();
				mesh->DeepCopy(meshReader->GetOutput());
			
				delete meshReader;
			}
			
		
			// extracting cells
			vtkGenericCell* cell = vtkGenericCell::New();
			vtkIdList* idList;
			vtkTetra* tetra = vtkTetra::New();
			vtkMeshQuality* quality = vtkMeshQuality::New();
		
			int count = 0;
			_edgeRatio.set_size(mesh->GetNumberOfCells());
			
			for (int y=0; y<mesh->GetNumberOfCells(); y++){
				
				// get mesh element
				mesh->GetCell(y, cell);
				idList = cell->GetPointIds();
				
				// put element vertices in the tetra
				double pt0[3]; mesh->GetPoint(idList->GetId(0),pt0); tetra->GetPoints()->SetPoint(0, pt0);
				double pt1[3]; mesh->GetPoint(idList->GetId(1),pt1); tetra->GetPoints()->SetPoint(1, pt1);
				double pt2[3]; mesh->GetPoint(idList->GetId(2),pt2); tetra->GetPoints()->SetPoint(2, pt2);
				double pt3[3]; mesh->GetPoint(idList->GetId(3),pt3); tetra->GetPoints()->SetPoint(3, pt3);

				// put element connections in the tetra
				tetra->GetPointIds()->SetId(0,0);
				tetra->GetPointIds()->SetId(1,1);
				tetra->GetPointIds()->SetId(2,2);
				tetra->GetPointIds()->SetId(3,3);

				// calculate the edge ratio
				_edgeRatio(y) = quality->TetEdgeRatio(tetra);
				if (quality->TetEdgeRatio(tetra) <= 0){
					count ++;	
				}
			}	

			// edgeRatio file name
			QString fileName = _meshFileNames[i];
			if (fileType == true)
				fileName.replace(QString(".inp"),QString("_edgeRatio.txt"));
			else if (fileType == false)
				fileName.replace(QString(".cdb"),QString("_edgeRatio.txt"));
			// writing		
			VnlWriterVector* vectorWriter = new VnlWriterVector;
			vectorWriter->SetVnlVector(_edgeRatio);
			vectorWriter->SetFileName(fileName);
			vectorWriter->Update();
	
			// extracting mesh id name (original meshes)
			/*QString temp = _meshFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("\\")+1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("/")+1);
			}*/
			// extracting mesh id name (created instances)
			QString temp = _meshFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(0, temp.lastIndexOf("\\")+10);
				if (fileType == true)
					temp.remove(".inp");
				else if (fileType == false)
					temp.remove(".cdb");

			}
			else {
				temp.remove(0, temp.lastIndexOf("/")+10);
				if (fileType == true)
					temp.remove(".inp");
				else if (fileType == false)
					temp.remove(".cdb");
			}	
			
			// give it to edge ratio quality
			_edgeRatioQuality(0,i) = temp.toDouble();
			
			/*
			// average of the edge ratios
			std::cout << "lenght: " << _edgeRatio.size() << std::endl;
			std::cout << "edge ratio mean: " << _edgeRatio.mean() << std::endl;
			_edgeRatioQuality(1,i) = _edgeRatio.mean();

			// standard deviation of the edge ratio
			double stdDev = 0.0;
			for (int a=0; a<_edgeRatio.size(); a++){
				stdDev += (_edgeRatio(a) - _edgeRatio.mean())*(_edgeRatio(a) - _edgeRatio.mean());
			}
			stdDev /= (_edgeRatio.size()-1);
			stdDev = sqrt (stdDev);
			_edgeRatioQuality(2,i) = stdDev;
			std::cout << "edge ratio std dev: " << stdDev << std::endl;
			*/
			quality->Delete();
			mesh->Delete();
			delete vectorWriter;
			}
			
					
		// write edge ratio quality matrix
		if (folderName.lastIndexOf("/") == -1){
			folderName.remove(folderName.lastIndexOf("\\")+1,folderName.size()-1);
			}
		else {
			folderName.remove(folderName.lastIndexOf("/")+1,folderName.size()-1);
		}
		folderName.append("EdgeRatioQuality.txt");
		VnlWriterMatrix* writer = new VnlWriterMatrix; 
		writer->SetVnlMatrix(_edgeRatioQuality.transpose());
		writer->SetFileName(folderName);
		writer->MatrixShapeUpdate();
		delete writer;

		std::cout << "files written" << std::endl;
	}


	void MeshQuality::CalculateMinAngle(){

		// define if the meshes are from abaqus (.inp) or ansys (.cdb)
		bool fileType;
		if (_meshFileNames[0].endsWith(".inp"))
			fileType = true;
		else if (_meshFileNames[0].endsWith(".cdb"))
			fileType = false;
		else
			std::cout << "input files not supported" << std::endl;

		// folder path
		QString folderName = _meshFileNames[0];
		
		// min angle quality
		_minAngleQuality.set_size(3, _meshFileNames.size()); // mesh id, min angle average, min angle std dev


		// for each mesh
		for (int i=0; i<_meshFileNames.size(); i++){

			vtkPolyData* mesh = vtkPolyData::New();
			
			// if abaqus 
			if (fileType == true){
				// read mesh
				MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
				meshReader->SetFileName(_meshFileNames[i]);
				std::cout << _meshFileNames[i].toAscii().data() << std::endl; 
				meshReader->Update();
				mesh->DeepCopy(meshReader->GetOutput());
			}
			
			// if ansys
			else if (fileType == false){
				// read mesh
				MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
				meshReader->SetFileName(_meshFileNames[i]);
				std::cout << _meshFileNames[i].toAscii().data() << std::endl; 
				meshReader->Update();
				mesh->DeepCopy(meshReader->GetOutput());
			}
			
		
			// extracting cells
			vtkGenericCell* cell = vtkGenericCell::New();
			vtkIdList* idList;
			vtkTetra* tetra = vtkTetra::New();
			vtkMeshQuality* quality = vtkMeshQuality::New();
		
			int count = 0;
			_minAngle.set_size(mesh->GetNumberOfCells());
			
			for (int y=0; y<mesh->GetNumberOfCells(); y++){
				
				// get mesh element
				mesh->GetCell(y, cell);
				idList = cell->GetPointIds();
				
				// put element vertices in the tetra
				double pt0[3]; mesh->GetPoint(idList->GetId(0),pt0); tetra->GetPoints()->SetPoint(0, pt0);
				double pt1[3]; mesh->GetPoint(idList->GetId(1),pt1); tetra->GetPoints()->SetPoint(1, pt1);
				double pt2[3]; mesh->GetPoint(idList->GetId(2),pt2); tetra->GetPoints()->SetPoint(2, pt2);
				double pt3[3]; mesh->GetPoint(idList->GetId(3),pt3); tetra->GetPoints()->SetPoint(3, pt3);

				// put element connections in the tetra
				tetra->GetPointIds()->SetId(0,0);
				tetra->GetPointIds()->SetId(1,1);
				tetra->GetPointIds()->SetId(2,2);
				tetra->GetPointIds()->SetId(3,3);

				// calculate the min angle
				_minAngle(y) = quality->TetMinAngle(tetra);
				if (quality->TetMinAngle(tetra) <= 0){
					count ++;	
				}
			}	

			// min angle file name
			QString fileName = _meshFileNames[i];
			if (fileType == true)
				fileName.replace(QString(".inp"),QString("_minAngle.txt"));
			else if (fileType == false)
				fileName.replace(QString(".cdb"),QString("_minAngle.txt"));
			// writing		
			VnlWriterVector* vectorWriter = new VnlWriterVector;
			vectorWriter->SetVnlVector(_minAngle);
			vectorWriter->SetFileName(fileName);
			vectorWriter->Update();
	
			// extracting mesh id name (original meshes)
			/*QString temp = _meshFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("\\")+1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("/")+1);
			}*/
			// extracting mesh id name (created instances)
			QString temp = _meshFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(0, temp.lastIndexOf("\\")+10);
				if (fileType == true)
					temp.remove(".inp");
				else if (fileType == false)
					temp.remove(".cdb");

			}
			else {
				temp.remove(0, temp.lastIndexOf("/")+10);
				if (fileType == true)
					temp.remove(".inp");
				else if (fileType == false)
					temp.remove(".cdb");
			}	
			
			// give it to min angle quality
			_minAngleQuality(0,i) = temp.toDouble();
			
			/*
			// average of the min angles
			std::cout << "lenght: " << _minAngle.size() << std::endl;
			std::cout << "min angle mean: " << _minAngle.mean() << std::endl;
			_minAngleQuality(1,i) = _minAngle.mean();

			// standard deviation of the min angle
			double stdDev = 0.0;
			for (int a=0; a<_minAngle.size(); a++){
				stdDev += (_minAngle(a) - _minAngle.mean())*(_minAngle(a) - _minAngle.mean());
			}
			stdDev /= (_minAngle.size()-1);
			stdDev = sqrt (stdDev);
			_minAngleQuality(2,i) = stdDev;
			std::cout << "min angle std dev: " << stdDev << std::endl;
			*/

			quality->Delete();
			mesh->Delete();
			delete vectorWriter;
			}
					
		// write min angle quality matrix
		if (folderName.lastIndexOf("/") == -1){
			folderName.remove(folderName.lastIndexOf("\\")+1,folderName.size()-1);
			}
		else {
			folderName.remove(folderName.lastIndexOf("/")+1,folderName.size()-1);
		}
		folderName.append("MinAngleQuality.txt");
		VnlWriterMatrix* writer = new VnlWriterMatrix; 
		writer->SetVnlMatrix(_minAngleQuality.transpose());
		writer->SetFileName(folderName);
		writer->MatrixShapeUpdate();
		delete writer;

		std::cout << "files written" << std::endl;
	}



}