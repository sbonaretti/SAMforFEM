/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PCAMesh.h>

#include <time.h>
#include <vector>

#include <FemAssignerNodes.h>
#include <ImageHandler.h>
#include <MeshReaderAbaqus.h>
#include <MeshReaderMorpherVolume.h>
#include <MeshWriterAbaqus.h>
#include <MeshWriterAnsys.h>
#include <PointWriterXyz.h>
#include <VnlReaderEValues.h>
#include <VnlReaderMatrix.h>
#include <VnlReaderVector.h>
#include <VnlWriterEValues.h>
#include <VnlWriterMatrix.h>
#include <VnlWriterVector.h>
#include <StatisticsBasics.h>
#include <StatisticsDistanceCalculator.h>

#include <vtkPCAAnalysisFilter.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>

#include <QDir>
#include <QFile>
#include <QString>
#include <QTextStream>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>

using namespace image;
using namespace mesh;
using namespace points;
using namespace fem;
using namespace vnl;
using namespace statistics;


namespace pca{

	// constructor
	PCAMesh::PCAMesh(){

		_nOfInstances = 0;
	}

	// destructor
	PCAMesh::~PCAMesh(){
	}

	
	// overwritten virtual functions
	
	void PCAMesh::ShapePCA(){

		
		// mesh topology
		vtkPolyData* outputMesh = vtkPolyData::New();
		QString fileExtension;

		// PCA calculator
		vtkPCAAnalysisFilter* PCACalculator = vtkPCAAnalysisFilter::New();
		PCACalculator->SetNumberOfInputs(_meshFileNames.size());
		
		
	
	/**** input file type ****/
		
		// .cdb
		if (_meshFileNames[0].endsWith(".cdb")){ // Ansys quadratic volume meshes
			fileExtension = (".cdb");
			for (int i=0; i<_meshFileNames.size(); i++){
				MeshReaderMorpherVolume* reader = new MeshReaderMorpherVolume;
				std::cout << "mesh: " << i+1 << ' ' << _meshFileNames[i].ascii() << std::endl;
				reader->SetFileName(_meshFileNames[i]);
				reader->Update();
				PCACalculator->SetInput(i, reader->GetOutput());
				// _outputMesh is used to preserve the mesh topology
				if (i==0)
					outputMesh = reader->GetOutput();
				// cleaning up
				delete reader;
						
			}
		}

		// .inp
		else if (_meshFileNames[0].endsWith(".inp")){ // Abaqus quadratic volume meshes
			fileExtension = (".inp");
			for (int i=0; i<_meshFileNames.size(); i++){
				MeshReaderAbaqus* reader = new MeshReaderAbaqus;
				std::cout << "mesh: " << _meshFileNames[i].ascii() << std::endl;
				reader->SetFileName(_meshFileNames[i]);
				reader->Update();
				PCACalculator->SetInput(i, reader->GetOutput());
				// _outputMesh is used to preserve the mesh topology
				if (i==0)
					outputMesh = reader->GetOutput();
				// cleaning up
				delete reader;
			}
		}
		
		// .stl
		else if (_meshFileNames[0].endsWith(".stl")){ // Surface linear meshes
			fileExtension = (".stl");
			for (int i=0; i<_meshFileNames.size(); i++){
				vtkSTLReader* reader = vtkSTLReader::New();
				std::cout << "mesh: " << _meshFileNames[i].ascii() << std::endl;
				reader->SetFileName(_meshFileNames[i].ascii());
				reader->Update();
				PCACalculator->SetInput(i, reader->GetOutput());
				// _outputMesh is used to preserve the mesh topology
				if (i==0)
					outputMesh = reader->GetOutput();
				// cleaning up
				//reader->Delete();
			}
		}
		
		else
			std::cout << "input files not supported" << std::endl;
	
		std::cout << "loaded " << _meshFileNames.size() << " meshes" << std::endl;

		
		// write dataset coordinates as matrix(to be used for the combined PCA - each column is a mesh)
		vnl_matrix<double> dataset;
		dataset.set_size(_meshFileNames.size(), outputMesh->GetNumberOfPoints()*3);
		for (int i=0; i<_meshFileNames.size(); i++){ // each row is a mesn
			for (int j=0; j<PCACalculator->GetInput(i)->GetNumberOfPoints(); j++){
				double pt[3];
				PCACalculator->GetInput(i)->GetPoint(j,pt);
				dataset(i,j*3) = pt[0];
				dataset(i,j*3+1) = pt[1];
				dataset(i,j*3+2) = pt[2];
			}
		}
		dataset = dataset.transpose(); // each column is a mesh

		// filename
		if ( _meshFileNames[0].lastIndexOf("/") == -1){
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}		
		_meshFileNames[0].append(QString("mesh shape dataset coordinates.txt"));
		
		// writing
		VnlWriterMatrix* vnlWriterMatrix = new VnlWriterMatrix;
		vnlWriterMatrix->SetFileName(_meshFileNames[0]);
		vnlWriterMatrix->SetVnlMatrix(dataset);
		//vnlWriterMatrix->Update();
		vnlWriterMatrix->MatrixShapeUpdate();


	/**** calculate PCA ****/
		std::cout << "calculating pca" << std::endl;
		double start = clock();
		PCACalculator->Update();
		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "computation time of shape PCA on " << _meshFileNames.size() << " meshes: " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);


	/**** eigenvalues ****/
		std::cout << "write eigenvalues" << std::endl;
		
		vtkFloatArray* eValues = PCACalculator->GetEvals();
		// evalues in vnl
		vnl_vector<double> vnlEValues;
		vnlEValues.set_size(eValues->GetNumberOfTuples());
		for (int i=0; i<vnlEValues.size(); i++)
			vnlEValues(i) = eValues->GetValue(i);

		// normalized evalues in vnl
		float intEigenValuesSum=0.0;
		for(unsigned int i= 0; i< eValues->GetNumberOfTuples(); i++ )
			intEigenValuesSum += eValues->GetValue(i);
		vnl_vector<double> vnlNormalizedEvalues;
		vnlNormalizedEvalues.set_size(vnlEValues.size());
		for (int i=0; i<vnlEValues.size(); i++)
			vnlNormalizedEvalues(i) = vnlEValues(i) / intEigenValuesSum * 100;
		
		// filename
		if ( _meshFileNames[0].lastIndexOf("/") == -1){
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}		
		_meshFileNames[0].append(QString("mesh shape pca eigenvalues.txt"));
		
		// writing
		VnlWriterEValues* vnlWriterEValues = new VnlWriterEValues;
		vnlWriterEValues->SetFileName(_meshFileNames[0]);
		vnlWriterEValues->SetVnlEValues(vnlEValues);
		vnlWriterEValues->SetVnlNormalizedEValues(vnlNormalizedEvalues);
		vnlWriterEValues->Update();


		
	/**** eigenvectors ****/
		std::cout << "write eigenvectors" << std::endl;
		
		// write eigenvectors coordinates as matrix(to be used for the combined PCA - each column is an eigenvector) 
		vnl_matrix<double> eigenvectors;
		eigenvectors.set_size(_meshFileNames.size(), outputMesh->GetNumberOfPoints()*3);
		for (int i=0; i<_meshFileNames.size(); i++){ // each row is an eigenvector
			for (int j=0; j<PCACalculator->GetOutput(i)->GetNumberOfPoints(); j++){
				double pt[3];
				PCACalculator->GetOutput(i)->GetPoint(j,pt);
				eigenvectors(i,j*3) = pt[0];
				eigenvectors(i,j*3+1) = pt[1];
				eigenvectors(i,j*3+2) = pt[2];
			}
		}
		eigenvectors = eigenvectors.transpose(); // each column is an eigenvector
		
		// filename
		if ( _meshFileNames[0].lastIndexOf("/") == -1){
				_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}
		_meshFileNames[0].append(QString("mesh shape pca eigenvectors.txt"));
		
		// writing
		VnlWriterMatrix* vnlWriterEigenVector = new VnlWriterMatrix;
		vnlWriterEigenVector->SetFileName(_meshFileNames[0]);
		vnlWriterEigenVector->SetVnlMatrix(eigenvectors);
		//vnlWriterEigenVector->Update();
		vnlWriterEigenVector->MatrixShapeUpdate();

				
	/**** average ****/ 
		std::cout << "write average" << std::endl;
		vtkFloatArray* b = vtkFloatArray::New();
		b->InsertNextValue(0.0);		
		PCACalculator->GetParameterisedShape(b, outputMesh);
		
		// to vnl
		vnl_vector<double> mean;
		mean.set_size(outputMesh->GetNumberOfPoints()*3);
		for (int i=0; i<outputMesh->GetNumberOfPoints(); i++){
			double pt[3];
			outputMesh->GetPoint(i,pt);
			mean(i*3) = pt[0];
			mean(i*3+1) = pt[1];
			mean(i*3+2) = pt[2];
		}

		// filename
		if ( _meshFileNames[0].lastIndexOf("/") == -1){
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}
		_meshFileNames[0].append(QString("mesh shape pca average.txt"));
		
		// writing
		VnlWriterVector *vnlWriteVector = new VnlWriterVector;
		vnlWriteVector->SetFileName(_meshFileNames[0]);
		vnlWriteVector->SetVnlVector(mean);
		vnlWriteVector->Update();
		
		// save as mesh format
		if ( _meshFileNames[0].lastIndexOf("/") == -1){
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}
		_meshFileNames[0].append(QString("mesh shape pca average"));
		_meshFileNames[0].append(fileExtension);
		
		if (fileExtension == (".cdb")){
			MeshWriterAnsys* writer = new MeshWriterAnsys;
			writer->SetFileName(_meshFileNames[0].ascii());
			writer->MeshOn();
			writer->SetMesh(outputMesh);
			writer->Update();
			delete writer;
		}
		if (fileExtension == (".inp")){
			MeshWriterAbaqus* writer = new MeshWriterAbaqus;
			writer->SetFileName(_meshFileNames[0].ascii());
			writer->MeshOn();
			writer->SetMesh(outputMesh);
			writer->Update();
			delete writer;
		}
		if (fileExtension == (".stl")){
			vtkSTLWriter* writer = vtkSTLWriter::New();
			writer->SetFileName(_meshFileNames[0].ascii());
			writer->SetInput(outputMesh);
			writer->Update();
			writer->Delete();
		}

/*
///////////////////// comparison with my code /////////////////////
		std::cout << std::endl;
		std::cout << "comparison with my code" << std::endl;
		_shapeDatasetCoordinatesFileName = ("C:/0.Data/test data/mesh pca/vtk - my code comparison/mesh shape dataset coordinates.txt");
		_nOfInstances = 10;
		double nOfPoints = 191532; 
		vnl_matrix<double> dataMatrix;
		dataMatrix.set_size(_nOfInstances, nOfPoints*3); // horizontal matrix !!! nOfPoints*3
		QFile intensityDataset(_shapeDatasetCoordinatesFileName);
		intensityDataset.open(QIODevice::ReadOnly | QIODevice::Text);
		QTextStream readIntensityDataset(&intensityDataset);
		for (int j=0; j<_nOfInstances; j++){
			for (int i=0; i<nOfPoints*3; i++){
				readIntensityDataset >> dataMatrix(j,i);
			}
		}
		intensityDataset.close();
		
		std::cout << "calculating pca" << std::endl;
		//start = clock();
		
		// compute PCA
		SetMatrix(dataMatrix);
		
		PrincipalComponentAnalysis();

		end = clock();
		total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "computation time of intensity PCA on " << _meshFileNames.size() << " meshes: " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);

		// average
		std::cout << "average" << std::endl;
		vnl_vector<double> mean2;
		mean2.set_size(nOfPoints*3);
		mean2 = GetMean();
		for (int i=0; i<9; i++)
			std::cout << "mean2" << mean2(i) << std::endl;
		
		// eigenvalues
		std::cout << "writing eigenvalues" << std::endl; 
		vnl_vector<double> eVal;
		eVal.set_size(_nOfInstances);
		eVal = GetEValues();
		std::cout << eVal << std::endl;

		// eigenvalues
		std::cout << "normalized eigenvalues" << std::endl; 
		vnl_vector<double> normVal;
		normVal.set_size(_nOfInstances);
		normVal = GetNormalizedEValues();
		std::cout << normVal << std::endl;

		// eigenvectors
		std::cout << "writing eigenvectors" << std::endl; 
		vnl_matrix<double> eVect;
		eVect.set_size(nOfPoints*3,_nOfInstances);
		eVect = GetEVectors();
		
		 
		// eigenvectors - save points as vtkPolyData in .vtk files
		vtkPolyDataWriter* polyDataWriter = vtkPolyDataWriter::New();

		for(int i = 0; i < _nOfModes; ++i) {
			polyDataWriter->SetInput(PCACalculator->GetOutput(i));
			if ( _meshFileNames[0].lastIndexOf("/") == -1){
				_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
			}
			else {
				_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
			}
			_meshFileNames[0].append(QString("shapeEigenvector"));
			_meshFileNames[0].append(QString("%1").arg(i+1));
			_meshFileNames[0].append(".vtk");
			polyDataWriter->SetFileName(_meshFileNames[0].ascii());
			polyDataWriter->Update();
		}
		
		// cleaning up
		//outputMesh->Delete();
		//PCACalculator->Delete();
		//eValues->Delete();
		//polyDataWriter->Delete();
		//b->Delete();
				
		
		// number of modes used to create the instances
		_nOfModes = _PCAAnalysisFilter->GetModesRequiredFor(_modelVariability); 
		std::cout << "number of modes for the " << _modelVariability*100 << " of variation: " << _nOfModes << std::endl;
			
		// instances
		if (_nOfInstances !=0){
			std::cout << "number of istances to create: " << _nOfInstances << std::endl;
			
			for (int i=0; i<_nOfInstances; i++){
				// create the paramenter array
				vtkFloatArray* b = vtkFloatArray::New();
				for (int a=0; a<_nOfModes; a++){
					b->InsertNextValue(((rand() %100) /100.0) *6 -3);
				}
			
				// create instance
				_PCAAnalysisFilter->GetParameterisedShape(b, _outputMesh);
				
				// save the instance
				if ( _meshFileNames[0].lastIndexOf("/") == -1){
				_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
				}
				else {
					_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
				}
				_meshFileNames[0].append(QString("instance"));
				for (int a=0; a<_nOfModes; a++){ 
					_meshFileNames[0].append(QString("%1").arg(b->GetValue(a)));
					if (a != _nOfModes-1)
						_meshFileNames[0].append(QString("_"));
				}
				_meshFileNames[0].append(fileExtension);
				std::cout << "creation of the instance " << i+1 << ": " << _meshFileNames[0].ascii() << std::endl;

				if (fileExtension == (".cdb")){
					MeshWriterAnsys* writer = new MeshWriterAnsys;
					writer->SetFileName(_meshFileNames[0].ascii());
					writer->SetInput(_outputMesh);
					writer->Update();
				}
				if (fileExtension == (".inp")){
					MeshWriterAbaqus* writer = new MeshWriterAbaqus;
					writer->SetFileName(_meshFileNames[0].ascii());
					writer->SetInput(_outputMesh);
					writer->Update();
				}
				if (fileExtension == (".stl")){
					vtkSTLWriter* writer = vtkSTLWriter::New();
					writer->SetFileName(_meshFileNames[0].ascii());
					writer->SetInput(_outputMesh);
					writer->Update();
				}

				
				// cleaning up
				b->Delete();
			}
		}
	*/
	}

	void PCAMesh::IntensityPCA(){

		std::cout << "!!! change image type for float to short in FemAssigner.h and uncomment greyLevelsExtractor->SetImage(imageHandler->GetImage());" << std::endl;
	
	/**** Get intensities from images ****/

		// data matrix of intensities
		vnl_matrix<double> intensities;
				
		for (int i=0; i<_meshFileNames.size(); i++){
			
			// read mesh (ansys .cdb)
			MeshReaderMorpherVolume* reader = new MeshReaderMorpherVolume;
			std::cout << std::endl;
			std::cout << "mesh " << i+1 << ": " << _meshFileNames[i].ascii() << std::endl;
			reader->SetFileName(_meshFileNames[i]);
			reader->Update();
			if (i==0){
				// data matrix size
				intensities.set_size(reader->GetOutput()->GetNumberOfPoints(), _meshFileNames.size());
			}
			// read correspondent image
			QString temp = _meshFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("\\")+1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("/")+1);
			}
			ImageHandler* imageHandler = new ImageHandler;
			for (int y=0; y<_imageFileNames.size(); y++){
				if (_imageFileNames[y].contains(temp)){
					std::cout << "image: " << _imageFileNames[y].ascii() << std::endl;
					imageHandler->SetImageFileName(_imageFileNames[y].ascii());
					imageHandler->MetafileReaderUpdate();
				}
			}
			// image extrusion
			std::cout << "extruding and getting grey levels" << std::endl;
			imageHandler->Extrusion();

			// extract nodes grey levels
			FemAssignerNodes* greyLevelsExtractor = new FemAssignerNodes;
			greyLevelsExtractor->SetMesh(reader->GetOutput());
//greyLevelsExtractor->SetImage(imageHandler->GetImage()); // uncomment it after changing the image type
			greyLevelsExtractor->GreyLevelAssignmentUpdate();
			vtkDoubleArray* greyLevels = vtkDoubleArray::New();
			greyLevels = greyLevelsExtractor->GetGreyLevels();

			// assign to the pca data matrix intensities
			for (int a=0; a<greyLevels->GetNumberOfTuples(); a++){
				intensities(a,i) = greyLevels->GetValue(a) * 1000; // grey level are divided by 1000 in the FemAssigner class
			}

			// cleaning up
			delete reader;
			delete imageHandler;
			delete greyLevelsExtractor;
		}

		// write dataset grey levels (to be used for the combined PCA - each column is the intensity of one mesh)
		std::cout << "write intensity dataset grey levels" << std::endl;
		if ( _meshFileNames[0].lastIndexOf("/") == -1){
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}		
		_meshFileNames[0].append(QString("mesh intensity dataset grey levels.txt"));

		VnlWriterMatrix* vnlWriterMatrix = new VnlWriterMatrix;
		vnlWriterMatrix->SetFileName(_meshFileNames[0]);
		vnlWriterMatrix->SetVnlMatrix(intensities);
		vnlWriterMatrix->MatrixShapeUpdate();


	/**** calculate PCA ****/
		std::cout << "calcule pca" << std::endl;
		double start = clock();
		
		SetMatrix(intensities);
		PrincipalComponentAnalysis();

		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "computation time of intensity PCA on " << _meshFileNames.size() << " meshes: " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);

		
	/**** eigenvalues ****/
		std::cout << "write eigenvalues" << std::endl;
		
		vnl_vector<double> eValues;
		eValues.set_size(intensities.rows());
		eValues = GetEValues();
		
		vnl_vector<double> normalizedEValues;
		normalizedEValues.set_size(intensities.rows());
		normalizedEValues = GetNormalizedEValues();

		if ( _meshFileNames[0].lastIndexOf("/") == -1){
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}		
		_meshFileNames[0].append(QString("mesh intensity pca eigenvalues.txt"));

		// writing
		VnlWriterEValues* vnlWriterEValues = new VnlWriterEValues;
		vnlWriterEValues->SetFileName(_meshFileNames[0]);
		vnlWriterEValues->SetVnlEValues(eValues);
		vnlWriterEValues->SetVnlNormalizedEValues(normalizedEValues);
		vnlWriterEValues->Update();
		
		
	/**** eigenvectors ****/
		std::cout << "write eigenvectors" << std::endl;
		
		// write eigenvectors coordinates as matrix(to be used for the combined PCA - each column is an eigenvector) 
		vnl_matrix<double> eVectors;
		int c = intensities.columns();
		int r = intensities.rows();
		eVectors.set_size(r,c);
		eVectors = GetEVectors();

		// filename
		if ( _meshFileNames[0].lastIndexOf("/") == -1){
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}		
		_meshFileNames[0].append(QString("mesh intensity pca eigenvectors.txt"));
		
		// writing
		VnlWriterMatrix* vnlWriterEigenVector = new VnlWriterMatrix;
		vnlWriterEigenVector->SetFileName(_meshFileNames[0]);
		vnlWriterEigenVector->SetVnlMatrix(eVectors);
		vnlWriterEigenVector->MatrixShapeUpdate();


	/**** average ****/ 
		std::cout << "write average" << std::endl; 
		vnl_vector<double> mean;
		mean.set_size(r);
		mean = GetMean();

		// filename
		if ( _meshFileNames[0].lastIndexOf("/") == -1){
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("\\")+1,_meshFileNames[0].size()-1);
		}
		else {
			_meshFileNames[0].remove(_meshFileNames[0].lastIndexOf("/")+1,_meshFileNames[0].size()-1);
		}		
		_meshFileNames[0].append(QString("mesh intensity pca average.txt"));
		
		// writing
		VnlWriterVector *vnlWriteVector = new VnlWriterVector;
		vnlWriteVector->SetFileName(_meshFileNames[0]);
		vnlWriteVector->SetVnlVector(mean);
		vnlWriteVector->Update();

	}

	void PCAMesh::CombinedPCA(){

			
	/**** load shape pca outputs ****/ 
		
		// dataset coordinates
 		std::cout << "shape dataset coordinates: " << _shapeDatasetCoordinatesFileName.toAscii().data() << std::endl;
		vnl_matrix<double> datasetCoordinates;
		VnlReaderMatrix* readerMatrix = new VnlReaderMatrix;
		readerMatrix->SetFileName(_shapeDatasetCoordinatesFileName);
		readerMatrix->Update();
		datasetCoordinates = readerMatrix->GetVnlMatrix();
		
		int nOfPoints = datasetCoordinates.rows();
		int nOfCoordinates = datasetCoordinates.rows() * 3.0;
		int nOfShapes = datasetCoordinates.columns();
			
		// average
		std::cout << "shape average: " << _shapeAverageFileName.toAscii().data() << std::endl;
		vnl_vector<double> shapeAverage;
		VnlReaderVector* readerVector = new VnlReaderVector;
		readerVector->SetFileName(_shapeAverageFileName);
		readerVector->Update();
		shapeAverage = readerVector->GetVnlVector();
		// create the average matrix to subtract
		vnl_matrix<double> shapeAverageMatrix;
		shapeAverageMatrix.set_size(datasetCoordinates.rows(), datasetCoordinates.cols());
		for (int a=0; a<shapeAverageMatrix.cols(); a++){
			shapeAverageMatrix.set_column(a, shapeAverage);
		}

		 
		// eigenvalues
		std::cout << "shape eigenvalues: " << _shapeEValuesFileName.toAscii().data() << std::endl;
		vnl_vector<double> shapeEigenValues;
		VnlReaderEValues* readerEValues = new VnlReaderEValues;
		readerEValues->SetFileName(_shapeEValuesFileName);
		readerEValues->Update();
		shapeEigenValues = readerEValues->GetVnlEValues();
			
		// eigenvectors
		std::cout << "shape eigenvectors: " << _shapeEVectorsFileName.toAscii().data() << std::endl;
		vnl_matrix<double> shapeEigenVectors;
		readerMatrix->SetFileName(_shapeEVectorsFileName);
		readerMatrix->Update();
		shapeEigenVectors = readerMatrix->GetVnlMatrix();
				

	/**** load intensity pca outputs ****/ 

		// dataset intensities
		std::cout << "intensity dataset coordinates: " << _intensityDatasetCoordinatesFileName.toAscii().data() << std::endl;
		vnl_matrix<double> datasetIntensities;
		readerMatrix->SetFileName(_intensityDatasetCoordinatesFileName);
		readerMatrix->Update();
		datasetIntensities = readerMatrix->GetVnlMatrix();
						
		// average
		std::cout << "intensity average: " << _intensityAverageFileName.toAscii().data() << std::endl;
		vnl_vector<double> intensityAverage;
		readerVector->SetFileName(_intensityAverageFileName);
		readerVector->Update();
		intensityAverage = readerVector->GetVnlVector();
		// create the average matrix to subtract
		vnl_matrix<double> intensityAverageMatrix;
		intensityAverageMatrix.set_size(datasetIntensities.rows(), datasetIntensities.cols());
		for (int a=0; a<intensityAverageMatrix.cols(); a++){
			intensityAverageMatrix.set_column(a, intensityAverage);
		}

		// eigenvalues
		std::cout << "intensity eigenvalues: " << _intensityEValuesFileName.toAscii().data() << std::endl;
		vnl_vector<double> intensityEigenValues;
		readerEValues->SetFileName(_intensityEValuesFileName);
		readerEValues->Update();
		intensityEigenValues = readerEValues->GetVnlEValues();
		
		// eigenvectors
		std::cout << "intensity eigenvectors: " << _intensityEVectorsFileName.toAscii().data() << std::endl;
		vnl_matrix<double> intensityEigenVectors;
		readerMatrix->SetFileName(_intensityEVectorsFileName);
		readerMatrix->Update();
		intensityEigenVectors = readerMatrix->GetVnlMatrix();
		

	/**** compute shape parameters ****/

		std::cout << "compute shape parameters" << std::endl;

		// creation of W
		double sumShape = 0.0;
		for (int i=0; i<shapeEigenValues.size(); i++){
			sumShape += shapeEigenValues(i);
		}

		double sumIntensity = 0.0;
		for (int i=0; i<intensityEigenValues.size(); i++){
			sumIntensity += intensityEigenValues(i);
		}
		
		double w = std::sqrt (sumIntensity / sumShape);///////////////
		std::cout << "w: " << w << std::endl;

		vnl_matrix<double> W;
		W.set_size(nOfShapes,nOfShapes);
		W.fill(0.0);
		W.fill_diagonal(w);
		

		// calculation of the shape parameter matrix
		vnl_matrix<double> bs;
		bs.set_size(nOfShapes, nOfShapes);

		bs = W *shapeEigenVectors.transpose() * (datasetCoordinates - shapeAverageMatrix);
		// pseudo-inverse instead of transpose
		//vnl_matrix<double> invertedShapeEigenVectors;
		//invertedShapeEigenVectors = vnl_matrix_inverse<double> ( shapeEigenVectors.transpose() * shapeEigenVectors) * shapeEigenVectors.transpose();
		//bs = W * invertedShapeEigenVectors * (datasetCoordinates - shapeAverageMatrix);
		
		
	/**** compute intensity parameters ****/
		std::cout << "compute intensity parameters" << std::endl;

		vnl_matrix<double> bg;
		bg.set_size(nOfShapes, nOfShapes);

		bg = intensityEigenVectors.transpose() * (datasetIntensities - intensityAverageMatrix);
		//pseudo-inverse instead of transpose
		//vnl_matrix<double> invertedIntensityEigenVectors;
		//invertedIntensityEigenVectors = vnl_matrix_inverse<double> ( intensityEigenVectors.transpose() * intensityEigenVectors) * intensityEigenVectors.transpose();
		//bg = invertedIntensityEigenVectors * (datasetIntensities - intensityAverageMatrix);
		

		
	/**** compute combined parameters ****/
		std::cout << "compute combined parameters" << std::endl;
		
		vnl_matrix<double> b;
		b.set_size(nOfShapes*2, nOfShapes);
		b.update(bs, 0,0);
		b.update(bg, nOfShapes,0);
		
		// writing b
		VnlWriterMatrix* writerMatrix = new VnlWriterMatrix;
		writerMatrix->SetFileName("b_meshes.txt");
		writerMatrix->SetVnlMatrix(b);
		writerMatrix->MatrixShapeUpdate();
		delete writerMatrix;

		

			
	/**** compute pca ****/
		std::cout << "compute combined pca" << std::endl;
		
		double start = clock();
				
		SetMatrix(b);
		PrincipalComponentAnalysis();

		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "computation time of intensity PCA on " << _meshFileNames.size() << " meshes: " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);

	/**** eigenvalues ****/
		std::cout << "write eigenvalues" << std::endl;

		vnl_vector<double> eValues;
		eValues.set_size(nOfShapes);
		eValues	= GetEValues();
		vnl_vector<double> normalizedEValues;
		normalizedEValues.set_size(nOfShapes);
		normalizedEValues = GetNormalizedEValues();

		std::cout << "eValues: " << eValues << std::endl;
		std::cout << "normalizedEValues: " << normalizedEValues << std::endl;
		
		// filename
		if ( _shapeDatasetCoordinatesFileName.lastIndexOf("/") == -1){
			_shapeDatasetCoordinatesFileName.remove(_shapeDatasetCoordinatesFileName.lastIndexOf("\\")+1,_shapeDatasetCoordinatesFileName.size()-1);
		}
		else {
			_shapeDatasetCoordinatesFileName.remove(_shapeDatasetCoordinatesFileName.lastIndexOf("/")+1,_shapeDatasetCoordinatesFileName.size()-1);
		}		
		_shapeDatasetCoordinatesFileName.append(QString("mesh combined pca eigenvalues.txt"));
		std::cout << _shapeDatasetCoordinatesFileName.toAscii().data() << std::endl;

		// writing
		VnlWriterEValues* vnlWriterEValues = new VnlWriterEValues;
		vnlWriterEValues->SetFileName(_shapeDatasetCoordinatesFileName);
		vnlWriterEValues->SetVnlEValues(eValues);
		vnlWriterEValues->SetVnlNormalizedEValues(normalizedEValues);
		vnlWriterEValues->Update();


	/**** eigenvectors ****/
		std::cout << "write eigenvectors" << std::endl;
		
		// write eigenvectors coordinates as matrix(to be used for the combined PCA - each column is an eigenvector) 
		vnl_matrix<double> eVectors;
		eVectors.set_size(nOfShapes*2,nOfShapes);
		eVectors = GetEVectors();

		// filename
		if ( _shapeDatasetCoordinatesFileName.lastIndexOf("/") == -1){
			_shapeDatasetCoordinatesFileName.remove(_shapeDatasetCoordinatesFileName.lastIndexOf("\\")+1,_shapeDatasetCoordinatesFileName.size()-1);
		}
		else {
			_shapeDatasetCoordinatesFileName.remove(_shapeDatasetCoordinatesFileName.lastIndexOf("/")+1,_shapeDatasetCoordinatesFileName.size()-1);
		}		
		_shapeDatasetCoordinatesFileName.append(QString("mesh combined pca eigenvectors.txt"));
		std::cout << _shapeDatasetCoordinatesFileName.toAscii().data() << std::endl;
		
		// writing
		VnlWriterMatrix* vnlWriterEigenVector = new VnlWriterMatrix;
		vnlWriterEigenVector->SetFileName(_shapeDatasetCoordinatesFileName);
		vnlWriterEigenVector->SetVnlMatrix(eVectors);
		vnlWriterEigenVector->MatrixShapeUpdate();


	/**** average ****/ 
		std::cout << "write average" << std::endl; 
		vnl_vector<double> mean;
		mean.set_size(nOfShapes*2);
		mean = GetMean();

		// filename
		if ( _shapeDatasetCoordinatesFileName.lastIndexOf("/") == -1){
			_shapeDatasetCoordinatesFileName.remove(_shapeDatasetCoordinatesFileName.lastIndexOf("\\")+1,_shapeDatasetCoordinatesFileName.size()-1);
		}
		else {
			_shapeDatasetCoordinatesFileName.remove(_shapeDatasetCoordinatesFileName.lastIndexOf("/")+1,_shapeDatasetCoordinatesFileName.size()-1);
		}		
		_shapeDatasetCoordinatesFileName.append(QString("mesh combined pca average.txt"));
		std::cout << _shapeDatasetCoordinatesFileName.toAscii().data() << std::endl;

		// writing
		VnlWriterVector *vnlWriteVector = new VnlWriterVector;
		vnlWriteVector->SetFileName(_shapeDatasetCoordinatesFileName);
		vnlWriteVector->SetVnlVector(mean);
		vnlWriteVector->Update();


		delete vnlWriterEValues;
		delete vnlWriterEigenVector;
		delete vnlWriteVector;


	}
	void PCAMesh::InstanceCreation(){

		// loading the weights of the new instances
		LoadWeightsForInstanceCreation();
		
		// loading from combined pca
		LoadWandNofMeshes();
		LoadCombinedEValues(); // to calculate the weight of the parameters
		LoadCombinedEVectors();
		
		// loading from shape pca
		LoadShapeAverage();
		LoadShapeEVectors();

		// loading from intensity pca
		LoadIntensityAverage();
		LoadIntensityEVectors();

		
	// INSTANCE CREATION
	
		for (int i=0; i<_nOfInstances; i++){
			std::cout << std::endl;
			std::cout << "instance number: " << i+1 << std::endl;
			double start = clock();
		
			// instance folder creation
			QDir dir; 
			QString temp = _outputFolder;
			QString boneNumber = QString("%1").arg(i+1);  
			temp.append("/bone ");
			temp.append(boneNumber);
			dir.mkdir(temp);


		/* PARAMETERS COMPUTATION */
			// create the combined parameters (both for shape and intensity)
			vnl_vector<double> c;
			c.set_size(_nOfModes);

			for (int j=0; j<_nOfModes; j++){
				c(j) = _weights(j,i) * std::sqrt(_combinedEValues(j));
				std::cout << "parameter " << j+1 << ": " << c(j)  << " = " << _weights(j,i) << " * " << std::sqrt(_combinedEValues(j)) << std::endl;
			}
		

		
		/* SHAPE COMPUTATION */
		
			// vnl shape
			std::cout << "computing the shape" << std::endl;
			vnl_vector<double> x;
			x.set_size(_shapeAverage.size());
			x = _shapeAverage + (_shapeEigenvectors * _Winv * _combinedEigenVectors.extract(_nOfOriginalMeshes,_nOfModes,0,0) * c);

			// vnl to polydata
			vtkPoints* points = vtkPoints::New();
			for (int a=0; a<x.size(); a=a+3 ){
				double pt[3];
				pt[0]=x(a); pt[1]=x(a+1); pt[2]=x(a+2); 
				points->InsertNextPoint(pt);
			}
			_shapeAverageMesh->SetPoints(points);
			


		/* INTENSITY COMPUTATION */
				
			// vnl intensity
			std::cout << "computing the intensity" << std::endl;
			vnl_vector<double> g;
			g.set_size(_intensityAverage.size());
			g = _intensityAverage + (_intensityEigenvectors * _combinedEigenVectors.extract(_nOfOriginalMeshes, _nOfModes, _nOfOriginalMeshes, 0) * c);


		
		/* WRITING NEW INSTANCE */
			std::cout << "writing files" << std::endl;

			// mesh
			QString instanceFileName = temp;
			instanceFileName.append("\\instance");
			instanceFileName.append ("_");
			instanceFileName.append(QString("%1").arg(i+1)); 
			instanceFileName.append(".cdb");
					
			MeshWriterAnsys* meshWriter = new MeshWriterAnsys;
			meshWriter->MeshOn();
			meshWriter->SetMesh(_shapeAverageMesh);
			meshWriter->SetFileName(instanceFileName);
			meshWriter->Update();
			
			// node coordinates
			instanceFileName.replace(QString(".cdb"), QString("_nodes.txt"));
			PointWriterXyz* pointWriter = new PointWriterXyz;
			pointWriter->SetInput(_shapeAverageMesh->GetPoints());
			pointWriter->SetFileName(instanceFileName);
			pointWriter->Update();
			
			// node intensities
			instanceFileName.replace(QString("_nodes.txt"), QString("_intensities.txt"));
			VnlWriterVector* vectorWriter = new VnlWriterVector;
			vectorWriter->SetFileName(instanceFileName);
			vectorWriter->SetVnlVector(g);
			vectorWriter->Update();
			
			std::cout << "instance " << i+1 << " created" << std::endl;
			double end = clock();
			double total = (end-start)/CLOCKS_PER_SEC;
			std::cout << "computation time for the creation of instance " << i+1 << " : " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);


			
			// delete
			delete meshWriter;
			delete vectorWriter;
		}
	}


	void PCAMesh::InstanceRecreation(){

		// load mode numbers to use
		LoadModeNumbers();
		std::cout << _modeNumbers << std::endl;
		
		// loading from combined pca
		LoadWandNofMeshes();
		LoadCombinedEValues(); // to calculate the weight of the parameters
		LoadCombinedEVectors();
		
		// loading from shape pca
		LoadShapeAverage();
		LoadShapeEVectors();

		// loading from intensity pca
		LoadIntensityAverage();
		LoadIntensityEVectors();

		// variables
		int nOfInstances = _instanceShapeFileNames.size();
		vnl_vector<double> parameterDistances;
		parameterDistances.set_size(nOfInstances);
		vnl_vector<double> shapeDistances;
		shapeDistances.set_size(nOfInstances);
		vnl_vector<double> intensityDistances;
		intensityDistances.set_size(nOfInstances);
		vnl_matrix<double> finalStatistics;
		finalStatistics.set_size(_modeNumbers.size(),7);

		
	// INSTANCE RECREATION
		std::cout << std::endl;

		for (int m=0; m<_modeNumbers.size(); m++){
			_nOfModes = _modeNumbers(m);
			std::cout << std::endl;
			std::cout << "number of modes: " << _nOfModes << std::endl;
			
			for (int i=0; i<nOfInstances; i++){
				std::cout << " " << std::endl;
				std::cout << "recreating instance number " << i+1 <<std::endl;
				
				double start = clock();
			
			/* PARAMETERS COMPUTATION */
				std::cout << "computing the parameters" << std::endl;

				// loading the shape to recreate
				LoadShapeToRecreate (_instanceShapeFileNames[i]);

				// shape parameters
				vnl_vector<double> bs;
				bs.set_size(_nOfOriginalMeshes);
				bs = _W *_shapeEigenvectors.transpose() * (_shapeToRecreate - _shapeAverage);


				// loading the intensity to recreate
				LoadIntensityToRecreate (_instanceShapeFileNames[i]);

				// intensity parameters
				vnl_vector<double> bg;
				bg.set_size(_nOfOriginalMeshes);
				bg = _intensityEigenvectors.transpose() * (_intensityToRecreate - _intensityAverage);

				// combined parameters
				vnl_vector<double> b;
				b.set_size(2*_nOfOriginalMeshes);
				b.update(bs,0);
				b.update(bg,_nOfOriginalMeshes);
				for (int j=0; j<_nOfModes; j++)
					std::cout << "b " << j+1 << ": " << b(j)  << " weight: " << b(j)/std::sqrt(_combinedEValues(j)) << std::endl;
						
				// cTilde
				vnl_vector<double> cTilde;
				cTilde = _combinedEigenVectors.transpose() * b;
							
				// bTilde
				vnl_vector<double> bTilde;
				bTilde = _combinedEigenVectors.extract(b.size(), _nOfModes, 0, 0) * cTilde.extract(_nOfModes, 0);

				
			/* SHAPE COMPUTATION */
				std::cout << "computing the shape" << std::endl;

				// bTildeShape
				vnl_vector<double> bTildeShape;
				bTildeShape = (1/_w) * bTilde.extract(_nOfModes, 0);

				// vnl shape
				vnl_vector<double> x;
				x = _shapeAverage + _shapeEigenvectors.extract(_shapeToRecreate.size(), _nOfModes, 0, 0) * bTildeShape;
				
				// vnl to vtkPolyData
				vtkPoints* points = vtkPoints::New();
				for (int a=0; a<x.size(); a=a+3 ){
					double pt[3];
					pt[0]=x(a); pt[1]=x(a+1); pt[2]=x(a+2); 
					points->InsertNextPoint(pt);
				}
				_meshToRecreate->SetPoints(points);

				// writing files
				WriteShapeRecreated(_instanceShapeFileNames[i]);



			/* INTENSITY COMPUTATION */
				std::cout << "computing the intensity" << std::endl;
							
				// bTildeIntensity
				vnl_vector<double> bTildeIntensity;
				bTildeIntensity = bTilde.extract(_nOfModes, _nOfOriginalMeshes);
				
				// vnl intensity
				vnl_vector<double> g;
				g = _intensityAverage + _intensityEigenvectors.extract(_intensityToRecreate.size(), _nOfModes, 0, 0) * bTildeIntensity;

				// writing file
				WriteIntensityRecreated(_instanceShapeFileNames[i], g);


				
			/* STATISTICS */
				StatisticsDistanceCalculator* distanceCalculator = new StatisticsDistanceCalculator;
				
				// parameters Euclidean distance
				distanceCalculator->SetVectorOne(b);
				distanceCalculator->SetVectorTwo(bTilde);
				distanceCalculator->CalculateEuclideanDistance();
				std::cout << "parameter distance: " << distanceCalculator->GetDistanceAverage() << std::endl;
				parameterDistances(i) = distanceCalculator->GetDistanceAverage();
				vnl_vector<double> parameterDistance = distanceCalculator->GetDistanceVector();
				QString temp = _instanceShapeFileNames[i];
				if (temp.lastIndexOf("/") == -1){
					temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
					temp.remove(0, temp.lastIndexOf("\\")+1);
				}
				else {
					temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
					temp.remove(0, temp.lastIndexOf("/")+1);
				}
				QString outputFileName = _outputFolder;
				outputFileName.append("\\");
				outputFileName.append(temp);
				outputFileName.append("_parameterDistances_mode_");
				outputFileName.append(QString("%1").arg(_nOfModes)); 
				outputFileName.append(".txt");
				VnlWriterVector* writerVector = new VnlWriterVector;
				writerVector->SetFileName(outputFileName);
				writerVector->SetVnlVector(parameterDistance);
				writerVector->Update();

				
				// shape	
				distanceCalculator->SetVectorOne(_shapeToRecreate);
				distanceCalculator->SetVectorTwo(x);
				distanceCalculator->CalculateEuclideanDistancePoints();
				std::cout << "shape distance: " << distanceCalculator->GetDistanceAverage() << std::endl;
				shapeDistances(i) = distanceCalculator->GetDistanceAverage();
				vnl_vector<double> shapeDistance = distanceCalculator->GetDistanceVector();
				outputFileName.replace(QString("_parameterDistances_mode_"), QString("_shapeDistances_mode_"));
				writerVector->SetFileName(outputFileName);
				writerVector->SetVnlVector(shapeDistance);
				writerVector->Update();
				
				// intensity	
				distanceCalculator->SetVectorOne(_intensityToRecreate);
				distanceCalculator->SetVectorTwo(g);
				distanceCalculator->CalculateEuclideanDistance();
				std::cout << "intensity distance: " << distanceCalculator->GetDistanceAverage() << std::endl;
				intensityDistances(i) = distanceCalculator->GetDistanceAverage();
				vnl_vector<double> intensityDistance = distanceCalculator->GetDistanceVector();
				outputFileName.replace(QString("_shapeDistances_mode_"), QString("_intensityDistances_mode_"));
				writerVector->SetFileName(outputFileName);
				writerVector->SetVnlVector(intensityDistance);
				writerVector->Update();

				double end = clock();
				double total = (end-start)/CLOCKS_PER_SEC;
				std::cout << "computation time for the recreation of the instance: " << total << " sec. about " << total/60 << " minutes" << std::endl;	


				// cleaning
				delete writerVector;
			}
			

		// FINAL STATISTICS
			std::cout << std::endl;
			std::cout << "final statistics" << std::endl; 	
			StatisticsBasics* stddevCalculator = new StatisticsBasics;
			
			// parameters
			std::cout << "parameters statistics:" << std::endl;
			std::cout << parameterDistances << std::endl;
			stddevCalculator->SetVector(parameterDistances);
			stddevCalculator->CalculateStandardDeviation();
			double parameterStdDev = stddevCalculator->GetStandardDeviation();
			std::cout << "distance: " << parameterDistances.mean() << " +- " << parameterStdDev << std::endl;
			
			// shape
			std::cout << "shape statistics:" << std::endl;
			std::cout << shapeDistances << std::endl;
			stddevCalculator->SetVector(shapeDistances);
			stddevCalculator->CalculateStandardDeviation();
			double shapeStdDev = stddevCalculator->GetStandardDeviation();
			std::cout << "distance: " << shapeDistances.mean() << " +- " << shapeStdDev << std::endl;

			// intensity
			std::cout << "intensity statistics:" << std::endl;
			std::cout << intensityDistances << std::endl;
			stddevCalculator->SetVector(intensityDistances);
			stddevCalculator->CalculateStandardDeviation();
			double intensityStdDev = stddevCalculator->GetStandardDeviation();
			std::cout << "distance: " << intensityDistances.mean() << " +- " << intensityStdDev << std::endl;
			
			// mode statistics
			finalStatistics(m,0) = _modeNumbers(m);
			finalStatistics(m,1) = parameterDistances.mean();
			finalStatistics(m,2) = parameterStdDev;
			finalStatistics(m,3) = shapeDistances.mean();
			finalStatistics(m,4) = shapeStdDev;
			finalStatistics(m,5) = intensityDistances.mean();
			finalStatistics(m,6) = intensityStdDev;
		}
	
		// write final statistics
		std::cout << "writing final statistics" << std::endl;
		QString outputFileName = _outputFolder;
		outputFileName.append("\\");
		outputFileName.append("recreation_statistics.txt");
		VnlWriterMatrix* writer = new VnlWriterMatrix;
		writer->SetFileName(outputFileName);
		writer->SetVnlMatrix(finalStatistics);
		writer->MatrixShapeUpdate();

		std::cout << "computations done" << std::endl;


		// cleaning
		delete writer;
	
	}


	void PCAMesh::TestMatrixPCA(){


		// read matrix
		VnlReaderMatrix* matrixReader = new VnlReaderMatrix;
		matrixReader->SetFileName(_shapeDatasetCoordinatesFileName);
		std::cout << _shapeDatasetCoordinatesFileName.toAscii().data() << std::endl;
		matrixReader->Update();
		
		/**** calculate PCA ****/
		std::cout << "calculating pca" << std::endl;
				
		SetMatrix(matrixReader->GetVnlMatrix());
		PrincipalComponentAnalysis();

				
	/**** eigenvalues ****/
		std::cout << "write eigenvalues" << std::endl;
		
		vnl_vector<double> eValues;
		eValues = GetEValues();
		
		vnl_vector<double> normalizedEValues;
		normalizedEValues = GetNormalizedEValues();

		// writing
		VnlWriterEValues* vnlWriterEValues = new VnlWriterEValues;
		vnlWriterEValues->SetFileName("eigenvalues.txt");
		vnlWriterEValues->SetVnlEValues(eValues);
		vnlWriterEValues->SetVnlNormalizedEValues(normalizedEValues);
		vnlWriterEValues->Update();
		
		
	/**** eigenvectors ****/
		std::cout << "write eigenvectors" << std::endl;
		
		// write eigenvectors coordinates as matrix(to be used for the combined PCA - each column is an eigenvector) 
		vnl_matrix<double> eVectors;
		eVectors = GetEVectors();

		// writing
		VnlWriterMatrix* vnlWriterEigenVector = new VnlWriterMatrix;
		vnlWriterEigenVector->SetFileName("eigenvectors.txt");
		vnlWriterEigenVector->SetVnlMatrix(eVectors);
		vnlWriterEigenVector->MatrixShapeUpdate();


	/**** average ****/ 
		std::cout << "write average" << std::endl; 
		vnl_vector<double> mean;
		mean = GetMean();

		// writing
		VnlWriterVector *vnlWriteVector = new VnlWriterVector;
		vnlWriteVector->SetFileName("average.txt");
		vnlWriteVector->SetVnlVector(mean);
		vnlWriteVector->Update();
	
	
	}
	void PCAMesh::ShapeInstanceRecreation(){

		int nOfInstances = _instanceShapeFileNames.size();

		// FROM SHAPE PCA //
		
		// load average mesh
		std::cout << "loading the average shape" << std::endl;
		std::cout << _shapeAverageFileName.toAscii().data() << std::endl;
		MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
		meshReader->SetFileName(_shapeAverageFileName);
		meshReader->Update();
		// to vnl
		vnl_vector<double> shapeAverage;
		shapeAverage.set_size(meshReader->GetOutput()->GetNumberOfPoints()*3);
		for (int i=0; i<meshReader->GetOutput()->GetNumberOfPoints(); i++){
			double pt[3];
			meshReader->GetOutput()->GetPoint(i,pt);
			shapeAverage(i*3) = pt[0];
			shapeAverage(i*3+1) = pt[1];
			shapeAverage(i*3+2) = pt[2];
		}

		//vectorReader->SetFileName(_shapeAverageFileName);//
		//vectorReader->Update();//
		//vnl_vector<double> shapeAverage;//
		//shapeAverage = vectorReader->GetVnlVector();//
									
		// load shape eigenvector 
		std::cout << "loading the shape eigenvectors" << std::endl;
		std::cout << _shapeEVectorsFileName.toAscii().data() << std::endl;
		VnlReaderMatrix* matrixReader = new VnlReaderMatrix; 
		matrixReader->SetFileName(_shapeEVectorsFileName);
		matrixReader->Update();
		vnl_matrix<double> shapeEigenvectors;
		shapeEigenvectors = matrixReader->GetVnlMatrix();

		// INSTANCES CALCULATION

		for (int i=0; i<nOfInstances; i++){
			std::cout << " " << std::endl;
					
			// load shape to recreate
			std::cout << "loading the shape to recreate" << std::endl;
			std::cout << _instanceShapeFileNames[i].toAscii().data() << std::endl;
			meshReader->SetFileName(_instanceShapeFileNames[i]);
			meshReader->Update();
			vtkPolyData* meshToRecreate= vtkPolyData::New();
			meshToRecreate = meshReader->GetOutput();
			// to vnl
			vnl_vector<double> shapeToRecreate;
			shapeToRecreate.set_size(meshReader->GetOutput()->GetNumberOfPoints()*3);
			for (int a=0; a<meshReader->GetOutput()->GetNumberOfPoints(); a++){
				double pt[3];
				meshReader->GetOutput()->GetPoint(a,pt);
				shapeToRecreate(a*3) = pt[0];
				shapeToRecreate(a*3+1) = pt[1];
				shapeToRecreate(a*3+2) = pt[2];
			}

			//vectorReader->SetFileName(_instanceShapeFileNames[i]);//
			//vectorReader->Update();//
			//vnl_vector<double> shapeToRecreate;//
			//shapeToRecreate = vectorReader->GetVnlVector();//

			// parameter calculation
			vnl_vector<double> b;
			b = shapeEigenvectors.transpose() * (shapeToRecreate-shapeAverage);


			// new shape calculation
			vnl_vector<double> x;
			x = shapeAverage + (shapeEigenvectors.extract(shapeAverage.size(),_nOfModes, 0,0) * b.extract(_nOfModes, 0));
 
			// vnl to vtkPolyData
			vtkPoints* points = vtkPoints::New();
			for (int a=0; a<x.size(); a=a+3 ){
				double pt[3];
				pt[0]=x(a); pt[1]=x(a+1); pt[2]=x(a+2); 
				points->InsertNextPoint(pt);
			}
			meshToRecreate->SetPoints(points);

			// WRITING FILES //
			std::cout << "writing files" << std::endl;

			// filename
			if (_instanceShapeFileNames[i].lastIndexOf("/") == -1){
				_instanceShapeFileNames[i].remove(_instanceShapeFileNames[i].lastIndexOf("\\")+6,_instanceShapeFileNames[i].size()-1);
				_instanceShapeFileNames[i].remove(0, _instanceShapeFileNames[i].lastIndexOf("\\")+1);
			}
			else {
				_instanceShapeFileNames[i].remove(_instanceShapeFileNames[i].lastIndexOf("/")+6,_instanceShapeFileNames[i].size()-1);
				_instanceShapeFileNames[i].remove(0, _instanceShapeFileNames[i].lastIndexOf("/")+1);
			}
			QString instanceFileName = _outputFolder;
			instanceFileName.append("\\");
			instanceFileName.append(_instanceShapeFileNames[i]);
			instanceFileName.append("_");
			instanceFileName.append(QString("%1").arg(_nOfModes)); //arg(i+1)
			instanceFileName.append("_modes.cdb");
			std::cout << instanceFileName.toAscii().data() << std::endl;
			
			// mesh
			MeshWriterAnsys* meshWriter = new MeshWriterAnsys;
			meshWriter->MeshOn();
			meshWriter->SetMesh(meshToRecreate);
			meshWriter->SetFileName(instanceFileName);
			meshWriter->Update();

			// original - new shape Euclidean distance
			StatisticsDistanceCalculator* distanceCalculator = new StatisticsDistanceCalculator;
			distanceCalculator->SetVectorOne(shapeToRecreate);
			distanceCalculator->SetVectorTwo(x);
			distanceCalculator->CalculateEuclideanDistance();
			std::cout << "distance average: " << distanceCalculator->GetDistanceAverage() << std::endl;
			std::cout << "distance standard deviation: " << distanceCalculator->GetDistanceStdDev() << std::endl;
			
		}
		
	}


	
	// private
	void PCAMesh::LoadShapeAverage(){

		std::cout << "loading the average shape" << std::endl;
		std::cout << _shapeAverageFileName.toAscii().data() << std::endl;
		
		// polydata
		MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
		meshReader->SetFileName(_shapeAverageFileName);
		meshReader->Update();
		_shapeAverageMesh = vtkPolyData::New();
		_shapeAverageMesh->DeepCopy(meshReader->GetOutput());

		// to vnl
		_shapeAverage.set_size(_shapeAverageMesh->GetNumberOfPoints()*3);
		for (int i=0; i<_shapeAverageMesh->GetNumberOfPoints(); i++){
			double pt[3];
			meshReader->GetOutput()->GetPoint(i,pt);
			_shapeAverage(i*3) = pt[0];
			_shapeAverage(i*3+1) = pt[1];
			_shapeAverage(i*3+2) = pt[2];
		}

		delete meshReader;

	
	}
	void PCAMesh::LoadShapeEValues(){}
	void PCAMesh::LoadShapeEVectors(){

		std::cout << "loading the shape eigenvectors" << std::endl;
		std::cout << _shapeEVectorsFileName.toAscii().data() << std::endl;
		VnlReaderMatrix* matrixReader = new VnlReaderMatrix; 
		matrixReader->SetFileName(_shapeEVectorsFileName);
		matrixReader->Update();
		_shapeEigenvectors = matrixReader->GetVnlMatrix();

		delete matrixReader;
	
	}
	void PCAMesh::LoadIntensityAverage(){
		
		std::cout << "loading the average intensity" << std::endl;
		std::cout << _intensityAverageFileName.toAscii().data() << std::endl;
		VnlReaderVector* vectorReader = new VnlReaderVector;
		vectorReader->SetFileName(_intensityAverageFileName);
		vectorReader->Update();
		_intensityAverage = vectorReader->GetVnlVector();

		delete vectorReader;
	
	}
	void PCAMesh::LoadIntensityEValues(){}
	void PCAMesh::LoadIntensityEVectors(){

		std::cout << "loading the intensity eigenvectors" << std::endl;
		std::cout << _intensityEVectorsFileName.toAscii().data() << std::endl;
		VnlReaderMatrix* matrixReader = new VnlReaderMatrix; 
		matrixReader->SetFileName(_intensityEVectorsFileName);
		matrixReader->Update();
		_intensityEigenvectors = matrixReader->GetVnlMatrix();

		delete matrixReader;
	}

	void PCAMesh::LoadCombinedEValues(){
	
		std::cout << "loading the combined eigenvalues" << std::endl;
		std::cout << _combinedEValuesFileName.toAscii().data() << std::endl;
		VnlReaderEValues* eValuesReader = new VnlReaderEValues;
		eValuesReader->SetFileName(_combinedEValuesFileName);
		eValuesReader->Update();
		_combinedEValues = eValuesReader->GetVnlEValues();

		delete eValuesReader;

	}
	void PCAMesh::LoadCombinedEVectors(){
	
		std::cout << "loading the combined eigenvectors" << std::endl;
		std::cout << _combinedEVectorsFileName.toAscii().data() << std::endl;
		VnlReaderMatrix* matrixReader = new VnlReaderMatrix; 
		matrixReader->SetFileName(_combinedEVectorsFileName);
		matrixReader->Update();
		_combinedEigenVectors = matrixReader->GetVnlMatrix();

		delete matrixReader;

	}
	void PCAMesh::LoadWandNofMeshes(){

		std::cout << "loading w and the number of original meshes" << std::endl;
		std::cout << _wFileName.toAscii().data() << std::endl;

		// number of original meshes
		VnlReaderVector* vectorReader = new VnlReaderVector;
		vectorReader->SetFileName(_wFileName);
		vectorReader->Update();
		vnl_vector<double> temp;
		temp = vectorReader->GetVnlVector();
		_nOfOriginalMeshes = temp(1);
		std::cout << "number of meshes in the original dataset: " << _nOfOriginalMeshes <<std::endl;
		
		// load w
		std::cout << "creating W" << std::endl;
		_w = temp(0);
		std::cout << "w: " << _w << std::endl;
		
		// create diagonal matrix W
		_W.set_size(_nOfOriginalMeshes);
		_W.fill(_w);
		
		// create diagonal matrix Winv
		_Winv.set_size(_nOfOriginalMeshes);
		_Winv.fill(1.0/_w);
		
		delete vectorReader;
			
	}

	void PCAMesh::LoadModeNumbers(){

		std::cout << "loading the mode numbers" << std::endl;
		std::cout << _modeNumbersFileName.toAscii().data() << std::endl;

		VnlReaderVector* reader = new VnlReaderVector;
		reader->SetFileName(_modeNumbersFileName);
		reader->Update();
		_modeNumbers = reader->GetVnlVector(); 

		delete reader;
	
	}

	void PCAMesh::LoadShapeToRecreate(QString fileName){

		std::cout << "loading the shape to recreate" << std::endl;
		std::cout << fileName.toAscii().data() << std::endl;
		
		// polydata
		MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
		meshReader->SetFileName(fileName);
		meshReader->Update();
		_meshToRecreate = vtkPolyData::New();
		_meshToRecreate->DeepCopy(meshReader->GetOutput());
		
		// to vnl
		_shapeToRecreate.set_size(_meshToRecreate->GetNumberOfPoints()*3);
		for (int a=0; a<_meshToRecreate->GetNumberOfPoints(); a++){
			double pt[3];
			_meshToRecreate->GetPoint(a,pt);
			_shapeToRecreate(a*3) = pt[0];
			_shapeToRecreate(a*3+1) = pt[1];
			_shapeToRecreate(a*3+2) = pt[2];
		}

		delete meshReader;
	
	}
	void PCAMesh::LoadIntensityToRecreate(QString fileName){
	
		std::cout << "loading the intensity to recreate" << std::endl;
		
		// loading the intensity corresponding to the shape
		VnlReaderVector* vectorReader = new VnlReaderVector;
		QString temp = fileName;
		if (temp.lastIndexOf("/") == -1){
			temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
			temp.remove(0, temp.lastIndexOf("\\")+1);
		}
		else {
			temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
			temp.remove(0, temp.lastIndexOf("/")+1);
		}
		for (int y=0; y<_instanceIntensityFileNames.size(); y++){
			if (_instanceIntensityFileNames[y].contains(temp)){
				vectorReader->SetFileName(_instanceIntensityFileNames[y]);
				std::cout << _instanceIntensityFileNames[y].ascii() << std::endl;
			}
		}

		vectorReader->Update();
		_intensityToRecreate = vectorReader->GetVnlVector();

		delete vectorReader;
	}
	void PCAMesh::LoadWeightsForInstanceCreation(){
	
		std::cout << "loading the combined weights" << std::endl;
		std::cout << _combinedWeightsFileName.toAscii().data() << std::endl;
		VnlReaderMatrix* matrixReader = new VnlReaderMatrix;
		matrixReader->SetFileName(_combinedWeightsFileName);
		matrixReader->Update();
		_weights.set_size(_nOfModes, _nOfInstances);
		std::cout << "number of modes: " << _weights.rows() << std::endl;
		std::cout << "number of instances to create: " << _weights.cols() << std::endl;
		_weights = matrixReader->GetVnlMatrix();

		delete matrixReader;
	
	}
	void PCAMesh::WriteShapeRecreated(QString fileName){

		std::cout << "writing the recreated shape" << std::endl;

		// filename
		if (fileName.lastIndexOf("/") == -1){
			fileName.remove(fileName.lastIndexOf("\\")+6,fileName.size()-1);
			fileName.remove(0, fileName.lastIndexOf("\\")+1);
		}
		else {
			fileName.remove(fileName.lastIndexOf("/")+6,fileName.size()-1);
			fileName.remove(0, fileName.lastIndexOf("/")+1);
		}
		
		QString instanceFileName = _outputFolder;
		instanceFileName.append("\\");
		instanceFileName.append(fileName);
		instanceFileName.append("_");
		instanceFileName.append(QString("%1").arg(_nOfModes)); //arg(i+1)
		instanceFileName.append("_modes.cdb");
		std::cout << instanceFileName.toAscii().data() << std::endl;
		
		// mesh
		MeshWriterAnsys* meshWriter = new MeshWriterAnsys;
		meshWriter->MeshOn();
		meshWriter->SetMesh(_meshToRecreate);
		meshWriter->SetFileName(instanceFileName);
		meshWriter->Update();

		// node coordinates
		instanceFileName.replace(QString(".cdb"), QString("_nodes.txt"));
		std::cout << instanceFileName.toAscii().data() << std::endl;
		PointWriterXyz* pointWriter = new PointWriterXyz;
		pointWriter->SetInput(_meshToRecreate->GetPoints());
		pointWriter->SetFileName(instanceFileName);
		pointWriter->Update();

		delete meshWriter;
		delete pointWriter;
	}

	void PCAMesh::WriteIntensityRecreated(QString fileName, vnl_vector<double> vector){

		std::cout << "writing the recreated intensity" << std::endl;

		// filename
		if (fileName.lastIndexOf("/") == -1){
			fileName.remove(fileName.lastIndexOf("\\")+6,fileName.size()-1);
			fileName.remove(0, fileName.lastIndexOf("\\")+1);
		}
		else {
			fileName.remove(fileName.lastIndexOf("/")+6,fileName.size()-1);
			fileName.remove(0, fileName.lastIndexOf("/")+1);
		}
		
		QString instanceFileName = _outputFolder;
		instanceFileName.append("\\");
		instanceFileName.append(fileName);
		instanceFileName.append("_");
		instanceFileName.append(QString("%1").arg(_nOfModes)); //arg(i+1)
		instanceFileName.append("_intensity.txt");
		std::cout << instanceFileName.toAscii().data() << std::endl;
		
		// write file
		VnlWriterVector* writer = new VnlWriterVector;
		writer->SetFileName(instanceFileName);
		writer->SetVnlVector(vector);
		writer->Update();

		delete writer;
	}

}

// PSEUDO-INVERSE CALCULATION //
/*	double start = clock();

// pseudo-inverse (Moore-Penrose method) for the shape parameters calculation (FiS * Ws^-1 * FiC,S)^-1
std::cout << "computing matrix inversion for shape parameters" << std::endl;
vnl_matrix<double> productS;
productS = shapeEigenvectors * Winv * combinedEigenVectors.extract(nOfOriginalMeshes, nOfOriginalMeshes, 0, 0);
vnl_matrix<double> invertedS;
invertedS = vnl_matrix_inverse<double> ( productS.transpose() * productS);
invertedS = invertedS * productS.transpose();

// pseudo-inverse (Moore-Penrose method) for the intensity parameters calculation (FiG * FiC,G)^-1
std::cout << "computing matrix inversion for intensity parameters" << std::endl;
vnl_matrix<double> productG;
productG = intensityEigenvectors * combinedEigenVectors.extract(nOfOriginalMeshes, nOfOriginalMeshes, nOfOriginalMeshes, 0);
vnl_matrix<double> invertedG;
invertedG = vnl_matrix_inverse<double> ( productG.transpose() * productG);
invertedG = invertedG * productG.transpose();

double end = clock();
double total = (end-start)/CLOCKS_PER_SEC;
std::cout << "computation time parameters " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);
		*/