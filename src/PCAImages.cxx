/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PCAImages.h>

#include <time.h>

#include <StatisticsDistanceCalculator.h>

#include <VnlReaderEValues.h>
#include <VnlReaderVector.h>
#include <VnlReaderMatrix.h>
#include <VnlWriterEValues.h>
#include <VnlWriterMatrix.h>
#include <VnlwriterVector.h>

#include <itkExponentialDeformationFieldImageFilter.h>
#include <itkFieldPCAShapeModelEstimator.h>
#include <itkImagePCAShapeModelEstimator.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkWarpImageFilter.h>

#include <QDir>

#include <vnl/vnl_vector.h>

#include <vtkPCAAnalysisFilter.h>//
#include <vtkPolyData.h>//
#include <vtkPoints.h>//


using namespace statistics;
using namespace vnl;


namespace pca{

	// constructor
	PCAImages::PCAImages(){
	}

	// destructor
	PCAImages::~PCAImages(){
	}

	
	// overwritten virtual functions
	void PCAImages::ShapePCA(){

		// PCA calculator
		typedef itk::FieldPCAShapeModelEstimator< FieldType,FieldType > FieldPCAShapeModelEstimatorType;
		FieldPCAShapeModelEstimatorType::Pointer PCACalculator = FieldPCAShapeModelEstimatorType::New();
		PCACalculator->SetNumberOfTrainingImages(_imageFileNames.size());
		PCACalculator->SetNumberOfPrincipalComponentsRequired(_imageFileNames.size());
		
		// load files in the calculator
		for(int i = 0; i < _imageFileNames.size(); ++i) {
			
			// load the vector image
			FieldReaderType::Pointer fieldReader = FieldReaderType::New();
			fieldReader->SetFileName(_imageFileNames[i].toAscii().data());
			std::cout << i+1 << ' ' << _imageFileNames[i].toAscii().data() << std::endl;
			fieldReader->Update();

			// put the vf in the calculator
			PCACalculator->SetInput(i, fieldReader->GetOutput());
		}
		std::cout << "loaded " << _imageFileNames.size() << " vector images" << std::endl;


	/**** calculate PCA ****/
		std::cout << "calculating PCA" << std::endl;
		double start = clock();
		PCACalculator->Update();
		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "computation time of PCA on " << _imageFileNames.size() << " vector images: " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);

	
	/**** eigenvalues ****/
		// evalues in vnl
		vnl_vector<double> vnlEValues;
		vnlEValues.set_size(PCACalculator->GetEigenValues().size());
		vnl_vector<double> vnlNormalizedEvalues;
		vnlNormalizedEvalues.set_size(PCACalculator->GetEigenValues().size());
		for (int i=0; i<vnlEValues.size(); i++){
			vnlEValues(i) = PCACalculator->GetEigenValues()[i];
			vnlNormalizedEvalues(i) = PCACalculator->GetEigenValuesNormalized()[i] * 100;
		}
		std::cout << "normalized eigenvalues:" << std::endl;
		std::cout << vnlNormalizedEvalues << std::endl;

		// filename
		std::cout << "write eigenvalues" << std::endl;
		if ( _imageFileNames[0].lastIndexOf("/") == -1){
			_imageFileNames[0].remove(_imageFileNames[0].lastIndexOf("\\")+1,_imageFileNames[0].size()-1);
		}
		else {
			_imageFileNames[0].remove(_imageFileNames[0].lastIndexOf("/")+1,_imageFileNames[0].size()-1);
		}		
		_imageFileNames[0].append(QString("shape model eigenvalues.txt"));
		std::cout << _imageFileNames[0].toAscii().data() << std::endl;
		
		// writing
		VnlWriterEValues* vnlWriterEValues = new VnlWriterEValues;
		vnlWriterEValues->SetFileName(_imageFileNames[0]);
		vnlWriterEValues->SetVnlEValues(vnlEValues);
		vnlWriterEValues->SetVnlNormalizedEValues(vnlNormalizedEvalues);
		vnlWriterEValues->Update();


	/**** eigenvectors ****/
		FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
		std::cout << "write eigenvectors" << std::endl;

		for(int i = 0; i < _imageFileNames.size(); ++i) {
			_imageFileNames[0].remove(_imageFileNames[0].lastIndexOf("\\")+1, _imageFileNames[0].size()-1);
			_imageFileNames[0].append("PC");
			_imageFileNames[0].append(QString("%1").arg(i+1));
			_imageFileNames[0].append("Shape.mhd");
			std::cout << _imageFileNames[0].toAscii().data() << std::endl;
			fieldWriter->SetFileName(_imageFileNames[0].toAscii().data());
			fieldWriter->SetInput(PCACalculator->GetOutput(i+1));
			fieldWriter->Update();
		}

		
	/**** average ****/
		std::cout << "write average" << std::endl;
		_imageFileNames[0].remove(_imageFileNames[0].lastIndexOf("\\")+1, _imageFileNames[0].size()-1);
		_imageFileNames[0].append("meanShape.mhd");
		std::cout << _imageFileNames[0].toAscii().data() << std::endl;
		fieldWriter->SetFileName(_imageFileNames[0].toAscii().data());
		fieldWriter->SetInput(PCACalculator->GetOutput(0));
		fieldWriter->Update();

		std::cout << "shape PCA computed" << std::endl;

		
		// cleaning
		delete vnlWriterEValues;

	
	}
	void PCAImages::IntensityPCA(){

		// PCA calculator
		typedef itk::ImagePCAShapeModelEstimator< ImageType,ImageType > ImagePCAShapeModelEstimatorType;
		ImagePCAShapeModelEstimatorType::Pointer PCACalculator = ImagePCAShapeModelEstimatorType::New();
		PCACalculator->SetNumberOfTrainingImages(_imageFileNames.size());
		PCACalculator->SetNumberOfPrincipalComponentsRequired(_imageFileNames.size());

		// load files in the calculator
		for(int i = 0; i < _imageFileNames.size(); ++i) {
			
			// load the image
			ImageReaderType::Pointer imageReader = ImageReaderType::New();
			imageReader->SetFileName(_imageFileNames[i].toAscii().data());
			std::cout << _imageFileNames[i].toAscii().data() << std::endl;
			imageReader->Update();
			
			// put the image in the calculator
			PCACalculator->SetInput(i, imageReader->GetOutput());
ImageRegionIteratorType c (imageReader->GetOutput(), imageReader->GetOutput()->GetRequestedRegion());
for (c.GoToBegin(); !c.IsAtEnd(); ++c)
	std::cout << c.Get() << std::endl; 
		}
		std::cout << "loaded " << _imageFileNames.size() << " images" << std::endl;

	
	/**** calculate PCA ****/
		std::cout << "calculating PCA" << std::endl;
		double start = clock();
		PCACalculator->Update();
		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "computation time of PCA on " << _imageFileNames.size() << " images: " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);


	/**** eigenvalues ****/
		// evalues sum
		float intEigenValuesSum=0.0;
		for(unsigned int i= 0; i< PCACalculator->GetEigenValues().size(); i++ )
    		intEigenValuesSum += PCACalculator->GetEigenValues()[i];
    	
		// evalues in vnl
		vnl_vector<double> eValues;
		eValues.set_size(PCACalculator->GetEigenValues().size());
		vnl_vector<double> normalizedEValues;
		normalizedEValues.set_size(PCACalculator->GetEigenValues().size());
		for(unsigned int i=0; i< PCACalculator->GetEigenValues().size(); i++ ){
			eValues(i) = PCACalculator->GetEigenValues()[i];
			normalizedEValues(i) = PCACalculator->GetEigenValues()[i]/intEigenValuesSum*100;
		}
		std::cout << "normalized eigenvalues:" << std::endl;
		std::cout << normalizedEValues << std::endl;
    	
		// filename
		std::cout << "write eigenvalues" << std::endl;	
		if ( _imageFileNames[0].lastIndexOf("/") == -1){
			_imageFileNames[0].remove(_imageFileNames[0].lastIndexOf("\\")+1,_imageFileNames[0].size()-1);
		}
		else {
			_imageFileNames[0].remove(_imageFileNames[0].lastIndexOf("/")+1,_imageFileNames[0].size()-1);
		}		
		_imageFileNames[0].append(QString("intensity model eigenvalues.txt"));
		std::cout << _imageFileNames[0].toAscii().data() << std::endl;

		// writing
		VnlWriterEValues* vnlWriterEValues = new VnlWriterEValues;
		vnlWriterEValues->SetFileName(_imageFileNames[0]);
		vnlWriterEValues->SetVnlEValues(eValues);
		vnlWriterEValues->SetVnlNormalizedEValues(normalizedEValues);
		vnlWriterEValues->Update();		

	
	/**** eigenvectors ****/
		typedef itk::ImageFileWriter< ImageType > ImageWriterType;
        ImageWriterType::Pointer imageWriter = ImageWriterType::New();
		std::cout << "write eigenvectors" << std::endl;

		for(int i = 0; i < _imageFileNames.size(); ++i) {
			_imageFileNames[0].remove(_imageFileNames[0].lastIndexOf("\\")+1, _imageFileNames[0].size()-1);
			_imageFileNames[0].append("PC");
			_imageFileNames[0].append(QString("%1").arg(i+1));
			_imageFileNames[0].append("Intensity.mhd");
			std::cout << _imageFileNames[0].toAscii().data() << std::endl;
			imageWriter->SetFileName(_imageFileNames[0].toAscii().data());
			imageWriter->SetInput(PCACalculator->GetOutput(i+1));
			imageWriter->Update();
ImageRegionIteratorType c (PCACalculator->GetOutput(i+1), PCACalculator->GetOutput(i+1)->GetRequestedRegion());
for (c.GoToBegin(); !c.IsAtEnd(); ++c)
	std::cout << c.Get() << std::endl;
		}


	/**** average ****/ 
		std::cout << "write average" << std::endl;
		_imageFileNames[0].remove(_imageFileNames[0].lastIndexOf("\\")+1, _imageFileNames[0].size()-1);
		_imageFileNames[0].append("meanIntensity.mhd");
		std::cout << _imageFileNames[0].toAscii().data() << std::endl;
		imageWriter->SetFileName(_imageFileNames[0].toAscii().data());
		imageWriter->SetInput(PCACalculator->GetOutput(0));
		imageWriter->Update();
ImageRegionIteratorType f (PCACalculator->GetOutput(0), PCACalculator->GetOutput(0)->GetRequestedRegion());
for (f.GoToBegin(); !f.IsAtEnd(); ++f)
	std::cout << f.Get() << std::endl;

		std::cout << "intensity PCA computed" << std::endl;

	}

	void PCAImages::CombinedPCA(){

	
	/**** LOADING FOR COMBINED PCA ****/
		// load w
		VnlReaderVector* wReader = new VnlReaderVector;
		wReader->SetFileName(_wFileName);
		wReader->Update();
		vnl_vector<double> temp;
		temp = wReader->GetVnlVector();
		double w = temp(0);
		std::cout << "combined weight: " << w <<std::endl;
		double nOfImages = temp(1);
		std::cout << "number of images in the original dataset: " << nOfImages <<std::endl;
	

	/**** LOADING FROM SHAPE PCA ****/
		// average shape
		FieldReaderType::Pointer fieldReaderMeanVelField  = FieldReaderType::New();
		fieldReaderMeanVelField->SetFileName(_shapeAverageFileName.toAscii().data());
		std::cout << "average shape: " << _shapeAverageFileName.toAscii().data() << std::endl;
		fieldReaderMeanVelField->Update();
	

	/**** LOADING FROM INTENSITY PCA ****/
		// average intensity 
		ImageReaderType::Pointer imageReaderMeanIntensity = ImageReaderType::New();
		imageReaderMeanIntensity->SetFileName(_intensityAverageFileName);
		std::cout << "average intensity: " << _intensityAverageFileName.toAscii().data() << std::endl;
		imageReaderMeanIntensity->Update();
		

		
	/**** COMBINED PCA ****/	
		// matrix of parameters
		vnl_matrix<double> b;
		b.set_size(nOfImages*2, nOfImages);
 		

		/** SHAPE PARAMETERS **/
		std::cout << "calculation of bs" << std::endl;
		
		for(int i = 0; i < nOfImages; ++i) {
			
			// read the i-th shape eigenvector
			FieldReaderType::Pointer eVectorReader  = FieldReaderType::New();
			QString temp = ("PC");
			temp.append(QString("%1").arg(i+1));
			temp.append("Shape.mhd");
			for (int y=0; y<_shapeEVectorsFileNames.size(); y++){
				if (_shapeEVectorsFileNames[y].contains(temp)){
					eVectorReader->SetFileName(_shapeEVectorsFileNames[y].toAscii().data());
					std::cout << _shapeEVectorsFileNames[y].toAscii().data() << std::endl;
				}
			}
			eVectorReader->Update();
			
			// multiply by w
			typedef itk::MultiplyByConstantImageFilter <FieldType, float, FieldType> MultiplyByConstantImageFilter;
			MultiplyByConstantImageFilter::Pointer multiplyByConstantImageFilter = MultiplyByConstantImageFilter::New();
			multiplyByConstantImageFilter->SetInput(eVectorReader->GetOutput());
			multiplyByConstantImageFilter->SetConstant(w);
			multiplyByConstantImageFilter->Update();
						
			// read the velocity fields 
			for(int j = 0; j < nOfImages; ++j) {
				// in the same order as the dataset for the calculation of the shape pca!!!!
				FieldReaderType::Pointer fieldReader  = FieldReaderType::New();
				std::cout << _velocityFieldFileNames[j].toAscii().data() << std::endl;
				fieldReader->SetFileName(_velocityFieldFileNames[j].toAscii().data());
				fieldReader->Update();
						
				// subtract meanVelField
				typedef itk::SubtractImageFilter< FieldType, FieldType, FieldType> SubtractFieldFilterType;
				SubtractFieldFilterType::Pointer subtractFieldFilter = SubtractFieldFilterType::New();
				subtractFieldFilter->SetInput1(fieldReader->GetOutput());
				subtractFieldFilter->SetInput2(fieldReaderMeanVelField->GetOutput());
				subtractFieldFilter->Update();
				
				// eigenvector * detrended instance (sum of all products)
				FieldRegionIteratorType it0 (multiplyByConstantImageFilter->GetOutput(), multiplyByConstantImageFilter->GetOutput()->GetRequestedRegion());
				FieldRegionIteratorType it1 (subtractFieldFilter->GetOutput(), subtractFieldFilter->GetOutput()->GetRequestedRegion());
				for (it0.GoToBegin(), it1.GoToBegin(); !it0.IsAtEnd(); ++it0,++it1) {
					it1.Set(it0.Get() * it1.Get()); // it1.Set() overwrites the detrended image, not the eigenvector! 
				}

				// sum of the three parameters
				double sum0 = 0.0; double sum1 = 0.0; double sum2 = 0.0;
				for (it1.GoToBegin(); !it1.IsAtEnd(); ++it1) {
					sum0 += it1.Get()[0];
					sum1 += it1.Get()[1];
					sum2 += it1.Get()[2];
				}

				//b(i,j) = sum0; b(i+nOfImages,j) = sum1; b(i+2*nOfImages,j) = sum2;
				b(i,j) = sum0; 
			}
			
		}


		
		/** INTENSITY PARAMETERS **/
		std::cout << "calculation of bg" << std::endl;
		
		for(int i = 0; i < nOfImages; ++i) {
			// read the i-th intensity eigenvector
			ImageReaderType::Pointer eVectorReader = ImageReaderType::New();
			QString temp = ("PC");
			temp.append(QString("%1").arg(i+1));
			temp.append("Intensity.mhd");
			for (int y=0; y<_intensityEVectorsFileNames.size(); y++){
				if (_intensityEVectorsFileNames[y].contains(temp)){
					eVectorReader->SetFileName(_intensityEVectorsFileNames[y]);
					std::cout << _intensityEVectorsFileNames[y].ascii() << std::endl;
				}
			}
			eVectorReader->Update();

			// read the images
			for(int j = 0; j < nOfImages; ++j) {
				// in the same order as the dataset for the calculation of the intensity pca!!!!
				ImageReaderType::Pointer imageReader = ImageReaderType::New();
				imageReader->SetFileName(_imageFileNames[j]);
				std::cout << _imageFileNames[j].toAscii().data() << std::endl;
				imageReader->Update();

				// subtract meanIntensity
				typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType> SubtractFieldFilterType;
				SubtractFieldFilterType::Pointer subtractFieldFilter = SubtractFieldFilterType::New();
				subtractFieldFilter->SetInput1(imageReader->GetOutput());
				subtractFieldFilter->SetInput2(imageReaderMeanIntensity->GetOutput());
				subtractFieldFilter->Update();
				
				// eigenvector * detrended instance (sum of all products)
				typedef itk::ImageRegionIterator< ImageType > ImageRegionIteratorType;
				ImageRegionIteratorType it0 (eVectorReader->GetOutput(), eVectorReader->GetOutput()->GetRequestedRegion());
				ImageRegionIteratorType it1 (subtractFieldFilter->GetOutput(), subtractFieldFilter->GetOutput()->GetRequestedRegion());
				for (it0.GoToBegin(), it1.GoToBegin(); !it0.IsAtEnd(); ++it0,++it1) {
					it1.Set(it0.Get() * it1.Get()); // it1.Set() overwrites the detrended image, not the eigenvector! 
				}

				// sum of the three parameters
				double sum = 0.0; 
				for (it1.GoToBegin(); !it1.IsAtEnd(); ++it1) {
					sum += it1.Get();
				}
				//b(i+3*nOfImages,j) = sum;
				b(i+nOfImages,j) = sum;
			}
		}
std::cout << "b" << std::endl;
std::cout << b << std::endl;
		
		// writing b
		VnlWriterMatrix* writerMatrix = new VnlWriterMatrix;
		writerMatrix->SetFileName("b_images.txt");
		writerMatrix->SetVnlMatrix(b);
		writerMatrix->MatrixShapeUpdate();
		delete writerMatrix;
		
		
		/** COMBINED PCA **/
		
		
		/**** compute pca ****/
		std::cout << "calculating combined PCA" << std::endl;
		double start = clock();
		SetMatrix(b);
		PrincipalComponentAnalysis();
		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "computation time of PCA: " << total << " sec. about " << total/60 << " minutes" << std::endl;	// fieldPCAShapeModelEstimator->Print(std::cout);

		
		/**** eigenvalues ****/
		vnl_vector<double> eValues;
		eValues.set_size(nOfImages);
		eValues	= GetEValues();
		vnl_vector<double> normalizedEValues;
		normalizedEValues.set_size(nOfImages);
		normalizedEValues = GetNormalizedEValues();
		std::cout << "normalized eigenvalues:" << std::endl;
		std::cout << normalizedEValues << std::endl;

		// filename
		std::cout << "write eigenvalues" << std::endl;
		QString eigenValuesFileName = _shapeAverageFileName;
		if ( eigenValuesFileName.lastIndexOf("/") == -1){
			eigenValuesFileName.remove(eigenValuesFileName.lastIndexOf("\\")+1,eigenValuesFileName.size()-1);
		}
		else {
			eigenValuesFileName.remove(eigenValuesFileName.lastIndexOf("/")+1,eigenValuesFileName.size()-1);
		}		
		eigenValuesFileName.append(QString("combined model eigenvalues.txt"));
		std::cout << eigenValuesFileName.toAscii().data() << std::endl;

		// writing
		VnlWriterEValues* vnlWriterEValues = new VnlWriterEValues;
		vnlWriterEValues->SetFileName(eigenValuesFileName);
		vnlWriterEValues->SetVnlEValues(eValues);
		vnlWriterEValues->SetVnlNormalizedEValues(normalizedEValues);
		vnlWriterEValues->Update();

		
		
		/**** eigenvectors ****/
		std::cout << "write eigenvectors" << std::endl;
		
		// write eigenvectors coordinates as matrix(to be used for the combined PCA - each column is an eigenvector) 
		vnl_matrix<double> eVectors;
		eVectors.set_size(nOfImages*2,nOfImages);
		eVectors = GetEVectors();

		// filename
		QString eigenVectorsFileName = _shapeAverageFileName;
		if ( eigenVectorsFileName.lastIndexOf("/") == -1){
			eigenVectorsFileName.remove(eigenVectorsFileName.lastIndexOf("\\")+1,eigenVectorsFileName.size()-1);
		}
		else {
			eigenVectorsFileName.remove(eigenVectorsFileName.lastIndexOf("/")+1,eigenVectorsFileName.size()-1);
		}		
		eigenVectorsFileName.append(QString("combined model eigenvectors.txt"));
		std::cout << eigenVectorsFileName.toAscii().data() << std::endl;
		
		// writing
		VnlWriterMatrix* vnlWriterEigenVector = new VnlWriterMatrix;
		vnlWriterEigenVector->SetFileName(eigenVectorsFileName);
		vnlWriterEigenVector->SetVnlMatrix(eVectors);
		vnlWriterEigenVector->MatrixShapeUpdate();


		/**** average ****/ 
		std::cout << "write average" << std::endl; 
		vnl_vector<double> mean;
		mean = GetMean();

		// filename
		QString averageFileName = _shapeAverageFileName;
		if ( averageFileName.lastIndexOf("/") == -1){
			averageFileName.remove(averageFileName.lastIndexOf("\\")+1,averageFileName.size()-1);
		}
		else {
			averageFileName.remove(averageFileName.lastIndexOf("/")+1,averageFileName.size()-1);
		}		
		averageFileName.append(QString("combined model average.txt"));
		std::cout << averageFileName.toAscii().data() << std::endl;

		// writing
		VnlWriterVector *vnlWriteVector = new VnlWriterVector;
		vnlWriteVector->SetFileName(averageFileName);
		vnlWriteVector->SetVnlVector(mean);
		vnlWriteVector->Update();


		/*
//////////////////
		vnl_matrix<double> bReduced;
		bReduced.set_size(6,3);
		//bReduced.update(b.extract(3,3,0,0),0,0);
		//bReduced.update(b.extract(3,3,9,0),3,0);
		bReduced.update(b, 0, 0); 

std::cout << "bReduced" << std::endl;
std::cout << bReduced << std::endl;
		SetMatrix(bReduced);
		PrincipalComponentAnalysis();

		vnl_matrix<double> reducedVectors;
		reducedVectors.set_size(6,3);
		reducedVectors = GetEVectors();
std::cout << "reducedVectors" << std::endl;
std::cout << reducedVectors << std::endl;

		vnl_matrix<double> cTilde;
		cTilde = reducedVectors.transpose() * bReduced;
		vnl_matrix<double> bTilde;
		bTilde = reducedVectors * cTilde;

std::cout << "bTilde" << std::endl;
std::cout << bTilde << std::endl;

		// comparison with vtk pca
		vtkPCAAnalysisFilter* PCACalculator = vtkPCAAnalysisFilter::New();
		PCACalculator->SetNumberOfInputs(3);
		
		double coord[3];
		
		std::cout << "first mesh" << std::endl;
		vtkPoints* pointsOne = vtkPoints::New();
		coord[0] = bReduced(0,0); coord[1] = bReduced(1,0); coord[2] = bReduced(2,0);
		std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
		pointsOne->InsertNextPoint(coord);
		coord[0] = bReduced(3,0); coord[1] = bReduced(4,0); coord[2] = bReduced(5,0);
		std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
		pointsOne->InsertNextPoint(coord);
		vtkPolyData* one = vtkPolyData::New();
		one->SetPoints(pointsOne);
		PCACalculator->SetInput(0, one);
				
		std::cout << "second mesh" << std::endl;
		vtkPoints* pointsTwo = vtkPoints::New();
		coord[0] = bReduced(0,1); coord[1] = bReduced(1,1); coord[2] = bReduced(2,1);
		std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
		pointsTwo->InsertNextPoint(coord);
		coord[0] = bReduced(3,1); coord[1] = bReduced(4,1); coord[2] = bReduced(5,1);
		std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
		pointsTwo->InsertNextPoint(coord);
		vtkPolyData* two = vtkPolyData::New();
		two->SetPoints(pointsTwo);
		PCACalculator->SetInput(1, two);
		
		std::cout << "third mesh" << std::endl;
		vtkPoints* pointsThree = vtkPoints::New();
		coord[0] = bReduced(0,2); coord[1] = bReduced(1,2); coord[2] = bReduced(2,2);
		std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
		pointsThree->InsertNextPoint(coord);
		coord[0] = bReduced(3,2); coord[1] = bReduced(4,2); coord[2] = bReduced(5,2);
		std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
		pointsThree->InsertNextPoint(coord);
		vtkPolyData* three = vtkPolyData::New();
		three->SetPoints(pointsThree);
		PCACalculator->SetInput(2, three);

		PCACalculator->Update();

		double pt[3];
		std::cout << "eigenvector1 " <<  std::endl;
		PCACalculator->GetOutput(0)->GetPoint(0,pt);
		std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
		PCACalculator->GetOutput(0)->GetPoint(1,pt);
		std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
		
		std::cout << "eigenvector2 " <<  std::endl;
		PCACalculator->GetOutput(1)->GetPoint(0,pt);
		std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
		PCACalculator->GetOutput(1)->GetPoint(1,pt);
		std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;

		std::cout << "eigenvector3 " <<  std::endl;
		PCACalculator->GetOutput(2)->GetPoint(0,pt);
		std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
		PCACalculator->GetOutput(2)->GetPoint(1,pt);
		std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;

		
		
		


///////////////////
*/
		// cleaning
		delete wReader;
		delete vnlWriterEValues;
		delete vnlWriterEigenVector;
		delete vnlWriteVector;
	}

	void PCAImages::InstanceCreation(){

		
	/**** LOADING FROM COMBINED PCA ****/
		// load w
		VnlReaderVector* wReader = new VnlReaderVector;
		wReader->SetFileName(_wFileName);
		wReader->Update();
		vnl_vector<double> temp;
		temp = wReader->GetVnlVector();
		double w = temp(0);
		std::cout << "w: " << w <<std::endl;
		std::cout << "1.0/w: " << 1.0/w << std::endl;
		double nOfOriginalImages = temp(1);
		std::cout << "number of images in the original dataset: " << nOfOriginalImages <<std::endl;
		
		// load combined eigenvalues
		std::cout << "loading the combined eigenvalues" << std::endl;
		VnlReaderEValues* eValuesReader = new VnlReaderEValues;
		eValuesReader->SetFileName(_combinedEValuesFileName);
		eValuesReader->Update();
		vnl_vector<double> combinedEValues;
		combinedEValues = eValuesReader->GetVnlEValues();
		
		// load combined weights for eigenvalues
		std::cout << "loading the combined weights" << std::endl;
		VnlReaderMatrix* matrixReader = new VnlReaderMatrix;
		matrixReader->SetFileName(_combinedWeightsFileName);
		matrixReader->Update();
		vnl_matrix<double> weights;
		weights.set_size(_nOfModes, _nOfInstances);
		std::cout << "number of modes: " << weights.rows() << std::endl;
		std::cout << "number of instances to create: " << weights.cols() << std::endl;
		weights = matrixReader->GetVnlMatrix();

		// load combined eigenvector matrix
		std::cout << "loading the combined eigenvectors" << std::endl;
		std::cout << _combinedEVectorsFileName.toAscii().data() << std::endl;
		matrixReader->SetFileName(_combinedEVectorsFileName);
		matrixReader->Update();
		vnl_matrix<double> combinedEigenVectors;
		combinedEigenVectors = matrixReader->GetVnlMatrix();
	

	/**** LOADING FROM SHAPE PCA ****/
		// load average velocity field
		std::cout << "loading the average shape" << std::endl;
		std::cout << _shapeAverageFileName.toAscii().data() << std::endl;
		FieldReaderType::Pointer averageVFreader = FieldReaderType::New();
		averageVFreader->SetFileName(_shapeAverageFileName.toAscii().data());
		averageVFreader->Update();


	/**** LOADING FROM INTENSITY PCA ****/
		// load average intensity
		std::cout << "loading the average intensity" << std::endl;
		std::cout << _intensityAverageFileName.toAscii().data() << std::endl;
		ImageReaderType::Pointer averageIntensityReader = ImageReaderType::New();
		averageIntensityReader->SetFileName(_intensityAverageFileName.toAscii().data());
		averageIntensityReader->Update();
		
	
	/**** INSTANCES CREATION ****/
	
		for (int i=0; i<_nOfInstances; i++){
		
			double start = clock(); 

			std::cout << std::endl;
			std::cout << "instance number " << i+1 << std::endl;

			// instance folder creation
			QDir dir; 
			QString temp = _outputFolder;
			QString boneNumber = QString("%1").arg(i+1);  
			temp.append("/bone ");
			temp.append(boneNumber);
			dir.mkdir(temp);
 
			// create the combined parameters (both for shape and intensity)
			vnl_vector<double> parameters;
			parameters.set_size(_nOfModes);
			for (int j=0; j<_nOfModes; j++){
				parameters(j) = weights(j,i) * std::sqrt(combinedEValues(j));
				std::cout << "parameters: " << j+1 << ": " << parameters(j)  << " = " << weights(j,i) << " * " << std::sqrt(combinedEValues(j)) << std::endl;
			}
	

		/**** SHAPE CREATION ****/
			std::cout << "Shape creation" << std::endl;
						
			// combined pca members
			vnl_vector<double> shapeParameters;
			shapeParameters.set_size(nOfOriginalImages*3);
			shapeParameters.fill(0.0);
			for (int a=0; a<_nOfModes; a++){
				for (int y=0; y<nOfOriginalImages; y++){
					shapeParameters(y) += (combinedEigenVectors(y,a) * parameters(a));
					//std::cout << shapeParameters(y) << ' ' << combinedEigenVectors(y,a) << ' ' << parameters(a) << std::endl;
				}
			}
			
			// vf to create
			FieldType::Pointer x = FieldType::New();
			x->SetRegions( averageVFreader->GetOutput()->GetRequestedRegion() );
			x->CopyInformation( averageVFreader->GetOutput() );
			x->Allocate();
			x->FillBuffer(0.0);				

			// for each shape eigenvector
			for (int a=0; a<nOfOriginalImages; a++){
						
				// load the shape eigenvector
				QString modeFileName = _shapeEVectorsFileNames[0];
				if (modeFileName.lastIndexOf("/") == -1)
					modeFileName.remove(modeFileName.lastIndexOf("\\")+1,modeFileName.size()-1);
				else 
					modeFileName.remove(modeFileName.lastIndexOf("/")+1,modeFileName.size()-1);
				modeFileName.append("PC");
				modeFileName.append(QString("%1").arg(a+1));
				modeFileName.append("Shape.mhd");
				std::cout << modeFileName.toAscii().data() << std::endl;
				
				FieldReaderType::Pointer shapeEVectorReader = FieldReaderType::New();
				shapeEVectorReader->SetFileName(modeFileName.toAscii().data());
				shapeEVectorReader->Update();
				
				// multiply the a-th shape eigenvector by 1/w
				std::cout << "multiply by 1/w: " << 1.0/w << std::endl;
				typedef itk::MultiplyByConstantImageFilter< FieldType, double, FieldType > MultiplyFilter;
				MultiplyFilter::Pointer multiplyFilter = MultiplyFilter::New();
				multiplyFilter->SetInput(shapeEVectorReader->GetOutput());
				multiplyFilter->SetConstant(1.0/w);
				multiplyFilter->Update();

				// multiply the result by the combined pca part
				std::cout << "multiply by the combined pca part" << std::endl;
				FieldRegionIteratorType it0 (multiplyFilter->GetOutput(), multiplyFilter->GetOutput()->GetRequestedRegion());
				FieldRegionIteratorType it1 (x, x->GetRequestedRegion());

				for (it0.GoToBegin(), it1.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1) {
					FieldType::PixelType temp;
					temp[0] = it1.Get()[0] + it0.Get()[0] * shapeParameters(a);
					temp[1] = it1.Get()[1] + it0.Get()[1] * shapeParameters(a);
					temp[2] = it1.Get()[2] + it0.Get()[2] * shapeParameters(a);
					it1.Set(temp);
					//std::cout << temp << std::endl;
				}
				
				
			}
			
			// sum up with the average velocity field
			std::cout << "sum the average" << std::endl;
			FieldRegionIteratorType it3 (x, x->GetRequestedRegion());
			FieldRegionIteratorType it4 (averageVFreader->GetOutput(), averageVFreader->GetOutput()->GetRequestedRegion() );
			for ( it3.GoToBegin(), it4.GoToBegin(); !it3.IsAtEnd(); ++it3, ++it4){
				it3.Set (it3.Get() + it4.Get());
				//std::cout << it3.Get() << std::endl;
			}

			// invert VF
			std::cout << "VF inversion" << std::endl;
			typedef itk::MultiplyByConstantImageFilter< FieldType, double, FieldType > MultiplyFilter;
			MultiplyFilter::Pointer multiplyFilter = MultiplyFilter::New();
			multiplyFilter->SetInput(x);
			multiplyFilter->SetConstant(-1);
			multiplyFilter->Update();
			
			// write SVF 
			std::cout << "write svf" << std::endl;
			QString temp4 = temp;
			temp4.append ("\\");
			temp4.append ("inverted_svf_");
			temp4.append(QString("%1").arg(i+1)); 
			//for (int n=0; n<_nOfModes; n++){
				//temp3.append ("_");
				//temp3.append(QString("%1").arg(weights(n,i))); 
			//}
			temp4.append (".mhd");
			typedef itk::ImageFileWriter< FieldType > FieldWriterType;
			FieldWriterType::Pointer fieldWriter2 = FieldWriterType::New();
			std::cout << temp4.toAscii().data() << std::endl;
			fieldWriter2->SetFileName(temp4.toAscii().data());
			fieldWriter2->SetInput(multiplyFilter->GetOutput());
			fieldWriter2->Update();
			
/*
			// VF to DVF 
			std::cout << "DVF creation (exponentiator)" << std::endl;
			typedef itk::ExponentialDeformationFieldImageFilter< FieldType,FieldType > ExponentialFieldFilterType;
			ExponentialFieldFilterType::Pointer exponentiator = ExponentialFieldFilterType::New();
			exponentiator->SetInput(multiplyFilter->GetOutput());
			exponentiator->Update();
			
			// write dvf 
			std::cout << "write dvf" << std::endl;
			QString temp3 = temp;
			temp3.append ("\\");
			temp3.append ("dvf_");
			temp3.append(QString("%1").arg(i+1)); 
			for (int n=0; n<_nOfModes; n++){
				//temp3.append ("_");
				//temp3.append(QString("%1").arg(weights(n,i))); 
			}
			temp3.append (".mhd");
			//typedef itk::ImageFileWriter< FieldType > FieldWriterType;
			FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
			std::cout << "temp3.toAscii().data()" << std::endl;
			fieldWriter->SetFileName(temp3.toAscii().data());
			fieldWriter->SetInput(exponentiator->GetOutput());
			fieldWriter->Update();
*/

	/**** INTENSITY CREATION ****/
/*			std::cout << "Intensity creation" << std::endl;
						
			// combined pca members
			vnl_vector<double> intensityParameters;
			intensityParameters.set_size(nOfOriginalImages);
			intensityParameters.fill(0.0);

			for (int a=0; a<_nOfModes; a++){
				for (int y=0; y<nOfOriginalImages; y++){
					intensityParameters(y) += (combinedEigenVectors((y + nOfOriginalImages),a) * parameters(a));
					//std::cout << intensityParameters(y) << ' ' << combinedEigenVectors((y + nOfOriginalImages),a) << ' ' << parameters(a) << std::endl;
				}
			}
						
			// image to create
			ImageType::Pointer g = ImageType::New();
			g->SetRegions( averageIntensityReader->GetOutput()->GetRequestedRegion() );
			g->CopyInformation( averageIntensityReader->GetOutput() );
			g->Allocate();
			g->FillBuffer(0.0);

			// for each intensity eigenvector
			for (int a=0; a<nOfOriginalImages; a++){

				// load the intensity eigenvector
				QString modeFileName = _intensityEVectorsFileNames[0];
				if (modeFileName.lastIndexOf("/") == -1)
					modeFileName.remove(modeFileName.lastIndexOf("\\")+1,modeFileName.size()-1);
				else 
					modeFileName.remove(modeFileName.lastIndexOf("/")+1,modeFileName.size()-1);
				modeFileName.append("PC");
				modeFileName.append(QString("%1").arg(a+1));
				modeFileName.append("Intensity.mhd");
				std::cout << modeFileName.toAscii().data() << std::endl;
				
				ImageReaderType::Pointer intensityEVectorReader = ImageReaderType::New();
				intensityEVectorReader->SetFileName(modeFileName.toAscii().data());
				intensityEVectorReader->Update();

				// multiply the result by the combined pca part
				std::cout << "multiply by the combined pca part" << std::endl;
				ImageRegionIteratorType it0 (intensityEVectorReader->GetOutput(), intensityEVectorReader->GetOutput()->GetRequestedRegion());
				ImageRegionIteratorType it1 (g, g->GetRequestedRegion());
				
				for (it0.GoToBegin(), it1.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1) {
					ImageType::PixelType temp;
					temp = it1.Get() + it0.Get() * intensityParameters(a);
					it1.Set(temp);
					//std::cout << temp << std::endl;
				}
			}
		
			// sum up with the average intensity
			std::cout << "sum the average" << std::endl;
			ImageRegionIteratorType it5 (g, g->GetRequestedRegion());
			ImageRegionIteratorType it6 (averageIntensityReader->GetOutput(), averageIntensityReader->GetOutput()->GetRequestedRegion() );
			for ( it5.GoToBegin(), it6.GoToBegin(); !it5.IsAtEnd(); ++it5, ++it6){
				it5.Set (it5.Get() + it6.Get());
				//std::cout << it5.Get() << std::endl;
			}
*/								
			
	/**** COMPOSING SHAPE AND INTENSITY ****/
/*			// warping intensity to shape
			typedef itk::WarpImageFilter< ImageType, ImageType, FieldType> WarpImageFilterType;
			WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();
			typedef WarpImageFilterType::CoordRepType CoordRepType;
			typedef itk::LinearInterpolateImageFunction< ImageType,CoordRepType > InterpolatorType;
			InterpolatorType::Pointer interpolator = InterpolatorType::New();
			
			warpImageFilter->SetInput(g);
			warpImageFilter->SetDeformationField(exponentiator->GetOutput());
			warpImageFilter->SetInterpolator(interpolator);
			warpImageFilter->SetOutputSpacing(g->GetSpacing());
			warpImageFilter->SetOutputOrigin(g->GetOrigin());
			warpImageFilter->SetEdgePaddingValue(-1024);
			warpImageFilter->Update();
			
			std::cout << "writing image warped to the reference" << std::endl;
			QString temp2 = temp;
			temp2.append ("\\");
			temp2.append ("intensityInReference_");
			temp2.append(QString("%1").arg(i+1)); 
			for (int n=0; n<_nOfModes; n++){
				//temp2.append ("_");
				//temp2.append(QString("%1").arg(weights(n,i))); 
			}
			temp2.append (".mhd");
			std::cout << temp2.toAscii().data() << std::endl;
			ImageWriterType::Pointer writer = ImageWriterType::New();
			writer->SetFileName(temp2.toAscii().data());
			writer->SetInput(g);
			writer->Update();

			std::cout << "writing image warped from the reference" << std::endl;
			temp.append ("\\");
			temp.append ("instance_");
			temp.append(QString("%1").arg(i+1)); 
			for (int n=0; n<_nOfModes; n++){
				//temp.append ("_");
				//temp.append(QString("%1").arg(weights(n,i))); 
			}
			temp.append (".mhd");
			std::cout << temp.toAscii().data() << std::endl;
			writer->SetFileName( temp.toAscii().data());
			writer->SetInput( warpImageFilter->GetOutput() );
			writer->Update();

			double end = clock();
			double total = (end-start)/CLOCKS_PER_SEC;
			std::cout << "Computation time for one instance: " << total << " sec. (" << total/60 << " min.)" << std::endl;
*/	
		}

		// cleaning
		delete wReader;
		delete eValuesReader;
		delete matrixReader;
		
	}


	void PCAImages::InstanceRecreation(){

	
		// loading the modes to use
		std::cout << "loading the mode numbers" << std::endl;
		std::cout << _modeNumbersFileName.toAscii().data() << std::endl;

		VnlReaderVector* reader = new VnlReaderVector;
		reader->SetFileName(_modeNumbersFileName);
		reader->Update();
		_modeNumbers = reader->GetVnlVector(); 
		delete reader;
		std::cout << _modeNumbers << std::endl;


		// allocating statistics  
		vnl_matrix<double> averageParameterDistances;
		averageParameterDistances.set_size(_instanceShapeFileNames.size(), _modeNumbers.size());
		vnl_matrix<double> averageShapeDistances;
		averageShapeDistances.set_size(_instanceShapeFileNames.size(), _modeNumbers.size());
		vnl_matrix<double> averageIntensityDistances;
		averageIntensityDistances.set_size(_instanceShapeFileNames.size(), _modeNumbers.size());
		

	/**** LOADING FROM COMBINED PCA ****/

		// load w
		std::cout << "creating W" << std::endl;
		std::cout << _wFileName.toAscii().data() << std::endl;
		VnlReaderVector* vectorReader = new VnlReaderVector;
		vectorReader->SetFileName(_wFileName);
		vectorReader->Update();
		vnl_vector<double> temp;
		temp = vectorReader->GetVnlVector();
		double w = temp(0);
		std::cout << "w: " << w << std::endl;
		// number of original meshes
		double nOfOriginalImages = temp(1);
		std::cout << "number of images in the original dataset: " << nOfOriginalImages <<std::endl;
		// number of instances to recreate
		int nOfInstances = _instanceShapeFileNames.size();
		std::cout << "number of images to recreate: " << nOfInstances <<std::endl;
		delete vectorReader;

		// load combined eigenvector matrix
		std::cout << "loading the combined eigenvectors" << std::endl;
		std::cout << _combinedEVectorsFileName.toAscii().data() << std::endl;
		VnlReaderMatrix* matrixReader = new VnlReaderMatrix; 
		matrixReader->SetFileName(_combinedEVectorsFileName);
		matrixReader->Update();
		vnl_matrix<double> combinedEigenVectors;
		combinedEigenVectors = matrixReader->GetVnlMatrix();
		delete matrixReader;
		
		// load combined eigenvalue vector (to calculate the weight of the parameters)
		std::cout << "loading the combined eigenvalues" << std::endl;
		std::cout << _combinedEValuesFileName.toAscii().data() << std::endl;
		VnlReaderEValues* eValuesReader = new VnlReaderEValues;
		eValuesReader->SetFileName(_combinedEValuesFileName);
		eValuesReader->Update();
		vnl_vector<double> combinedEValues;
		combinedEValues = eValuesReader->GetVnlEValues();
		delete eValuesReader;
			

	/**** LOADING FROM SHAPE PCA ****/

		// load average velocity field
		std::cout << "loading the average shape" << std::endl;
		FieldReaderType::Pointer averageVFreader = FieldReaderType::New();
		std::cout << _shapeAverageFileName.toAscii().data() << std::endl;
		averageVFreader->SetFileName(_shapeAverageFileName.toAscii().data());
		averageVFreader->Update();
//typedef itk::ImageRegionIterator< FieldType > FieldRegionIteratorType;
//FieldRegionIteratorType a (averageVFreader->GetOutput(), averageVFreader->GetOutput()->GetRequestedRegion());
//for (a.GoToBegin(); !a.IsAtEnd(); ++a)
//	std::cout << a.Get() << std::endl; 
		
			
	/**** LOADING FROM INTENSITY PCA ****/
		
		// load average intensity
		std::cout << "loading the average intensity" << std::endl;
		ImageReaderType::Pointer averageIntensityReader = ImageReaderType::New();
		std::cout << _intensityAverageFileName.toAscii().data() << std::endl;
		averageIntensityReader->SetFileName(_intensityAverageFileName.toAscii().data());
		averageIntensityReader->Update();




	/**** INSTANCE RECREATION ****/
	for (int m=0; m<_modeNumbers.size(); m++){
		_nOfModes = _modeNumbers(m);
		std::cout << std::endl;
		std::cout << "number of modes: " << _nOfModes << std::endl;
	
		
			for (int i=0; i<nOfInstances; i++){
			std::cout << std::endl;

			double start = clock();

			// load the velocity field to recreate
			std::cout << "loading the shape to recreate" << std::endl;
			FieldReaderType::Pointer fieldToRecreateReader = FieldReaderType::New();
			std::cout << _instanceShapeFileNames[i].toAscii().data() << std::endl;
			fieldToRecreateReader->SetFileName(_instanceShapeFileNames[i].toAscii().data());
			fieldToRecreateReader->Update();
	//FieldRegionIteratorType c (fieldToRecreateReader->GetOutput(), fieldToRecreateReader->GetOutput()->GetRequestedRegion());
	//for (c.GoToBegin(); !c.IsAtEnd(); ++c)
	//	std::cout << c.Get() << std::endl; 


			// load the intensity to recreate
			std::cout << "loading the intensity to recreate" << std::endl;
			ImageReaderType::Pointer imageToRecreateReader = ImageReaderType::New();
			QString temp = _instanceShapeFileNames[i];
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
					std::cout << _instanceIntensityFileNames[y].toAscii().data() << std::endl;
					imageToRecreateReader->SetFileName(_instanceIntensityFileNames[y].toAscii().data());
				}
			}
			imageToRecreateReader->Update();
			
			

		/**** COMPUTE PARAMETERS ****/
			
			/** shape parameters **/ 
			std::cout << "computing shape parameters" << std::endl;
			vnl_vector<double> bs;
			bs.set_size(nOfOriginalImages);
			
			// subtract meanVelField from instance
			typedef itk::SubtractImageFilter< FieldType, FieldType, FieldType> SubtractFieldFilterType;
			SubtractFieldFilterType::Pointer subtractFieldFilter = SubtractFieldFilterType::New();
			subtractFieldFilter->SetInput1(fieldToRecreateReader->GetOutput());
			subtractFieldFilter->SetInput2(averageVFreader->GetOutput());
			subtractFieldFilter->Update();
			
			for(int j = 0; j < nOfOriginalImages; ++j) {

				// empty image (sum of the weighted combined modes)
				FieldType::Pointer WBs = FieldType::New();
				WBs->SetRegions( averageVFreader->GetOutput()->GetRequestedRegion() );
				WBs->CopyInformation( averageVFreader->GetOutput() );
				WBs->Allocate();
				WBs->FillBuffer(0.0);
				
				// read the i-th shape eigenvector
				FieldReaderType::Pointer shapeEVectorReader = FieldReaderType::New();
				QString temp = ("PC");
				temp.append(QString("%1").arg(j+1));
				temp.append("Shape.mhd");
				for (int y=0; y<_shapeEVectorsFileNames.size(); y++){
					if (_shapeEVectorsFileNames[y].contains(temp)){
						shapeEVectorReader->SetFileName(_shapeEVectorsFileNames[y].toAscii().data());
						//std::cout << _shapeEVectorsFileNames[y].toAscii().data() << std::endl;
					}
				}
				shapeEVectorReader->Update();

	//FieldRegionIteratorType e (shapeEVectorReader->GetOutput(), shapeEVectorReader->GetOutput()->GetRequestedRegion());
	//for (e.GoToBegin(); !e.IsAtEnd(); ++e)
	//	std::cout << e.Get() << std::endl; 

				
				// multiply by w
				typedef itk::MultiplyByConstantImageFilter <FieldType, float, FieldType> MultiplyByConstantImageFilter;
				MultiplyByConstantImageFilter::Pointer multiplyByConstantImageFilter = MultiplyByConstantImageFilter::New();
				multiplyByConstantImageFilter->SetInput(shapeEVectorReader->GetOutput());
				multiplyByConstantImageFilter->SetConstant(w);
				multiplyByConstantImageFilter->Update();
				
				// eigenvector * detrended instance (sum of all products)
				FieldRegionIteratorType it0 (multiplyByConstantImageFilter->GetOutput(), multiplyByConstantImageFilter->GetOutput()->GetRequestedRegion());
				FieldRegionIteratorType it1 (subtractFieldFilter->GetOutput(), subtractFieldFilter->GetOutput()->GetRequestedRegion());
				FieldRegionIteratorType it2 (WBs, WBs->GetRequestedRegion());
				
				for (it0.GoToBegin(), it1.GoToBegin(), it2.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1, ++it2) {
					FieldType::PixelType value;
					value[0] = it0.Get()[0] * it1.Get()[0] + it2.Get()[0];
					value[1] = it0.Get()[1] * it1.Get()[1] + it2.Get()[1];
					value[2] = it0.Get()[2] * it1.Get()[2] + it2.Get()[2];
					it2.Set(value);
				}

				// sum of the three parameters
				double sum0 = 0.0;
				for (it2.GoToBegin(); !it2.IsAtEnd(); ++it2) {
					sum0 += it2.Get()[0];
					sum0 += it2.Get()[1];
					sum0 += it2.Get()[2];
				}
							
				bs(j) = sum0; //bs(j+nOfOriginalImages) = sum0; bs(j+2*nOfOriginalImages) = sum0;
			}		
	//std::cout << "bs" << std::endl;
	//std::cout << bs << std::endl;
					
			
			/** intensity parameters **/ 
			std::cout << "computing intensity parameters" << std::endl;
			vnl_vector<double> bg;
			bg.set_size(nOfOriginalImages);

			// subtract meanIntensity from instance
			typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType> SubtractImageFilterType;
			SubtractImageFilterType::Pointer subtractImageFilter = SubtractImageFilterType::New();
			subtractImageFilter->SetInput1(imageToRecreateReader->GetOutput());
			subtractImageFilter->SetInput2(averageIntensityReader->GetOutput());
			subtractImageFilter->Update();
			
			for(int j = 0; j < nOfOriginalImages; ++j) {
			
				// empty image (sum of the weighted combined modes)
				ImageType::Pointer Bg = ImageType::New();
				Bg->SetRegions( averageIntensityReader->GetOutput()->GetRequestedRegion() );
				Bg->CopyInformation( averageIntensityReader->GetOutput() );
				Bg->Allocate();
				Bg->FillBuffer(0.0);
				
				// read the i-th intensity eigenvector
				ImageReaderType::Pointer intensityEVectorReader = ImageReaderType::New();
				temp = ("PC");
				temp.append(QString("%1").arg(j+1));
				temp.append("Intensity.mhd");

				for (int y=0; y<_intensityEVectorsFileNames.size(); y++){
					if (_intensityEVectorsFileNames[y].contains(temp)){
						//std::cout << _intensityEVectorsFileNames[y].toAscii().data() << std::endl;
						intensityEVectorReader->SetFileName(_intensityEVectorsFileNames[y].toAscii().data());
					}
				}
				intensityEVectorReader->Update();
												
				// eigenvector * detrended instance (sum of all products)
				typedef itk::ImageRegionIterator< ImageType > RegionIteratorType;
				RegionIteratorType it3 (intensityEVectorReader->GetOutput(), intensityEVectorReader->GetOutput()->GetRequestedRegion());
				RegionIteratorType it4 (subtractImageFilter->GetOutput(), subtractImageFilter->GetOutput()->GetRequestedRegion());
				RegionIteratorType it5 (Bg, Bg->GetRequestedRegion());
				
				for (it3.GoToBegin(), it4.GoToBegin(), it5.GoToBegin(); !it3.IsAtEnd(); ++it3, ++it4, ++it5) {
					ImageType::PixelType value;
					value = it3.Get() * it4.Get() + it5.Get();
					it5.Set(value);
				}

				// sum of the three parameters
				double sum0 = 0.0;
				for (it5.GoToBegin(); !it5.IsAtEnd(); ++it5) {
					sum0 += it5.Get();
				}
							
				bg(j) = sum0; 
			}


			/** combined parameters **/ 
			std::cout << "computing combined parameters" << std::endl;
			vnl_vector<double> b;
			b.set_size(2*nOfOriginalImages);
			b.update(bs,0);
			b.update(bg,nOfOriginalImages);
			for (int j=0; j<_nOfModes; j++)
				std::cout << "b " << j+1 << ": " << b(j)  << " weight: " << b(j)/std::sqrt(combinedEValues(j)) << std::endl;
			/*
			vnl_matrix<double> combinedEig;
			combinedEig.set_size(nOfOriginalImages*2,nOfOriginalImages);
			combinedEig.update(combinedEigenVectors.extract(nOfOriginalImages,nOfOriginalImages,0,0),0,0);
			combinedEig.update(combinedEigenVectors.extract(nOfOriginalImages,nOfOriginalImages,nOfOriginalImages*3,0),nOfOriginalImages,0);
			combinedEigenVectors.set_size(nOfOriginalImages*2,nOfOriginalImages);
			combinedEigenVectors.update(combinedEig,0,0);
			*/
	//std::cout << "b" << std::endl;
	//std::cout << b << std::endl;
	//std::cout << "combinedEigenVectors" << std::endl;
	//std::cout << combinedEigenVectors << std::endl;

			// cTilde
			vnl_vector<double> cTilde;
			cTilde = combinedEigenVectors.transpose() * b;
	//std::cout << "cTilde" << std::endl;
	//std::cout << cTilde << std::endl;

			// bTilde
			vnl_vector<double> bTilde;
			bTilde = combinedEigenVectors.extract(b.size(), _nOfModes, 0, 0) * cTilde.extract(_nOfModes, 0);
	//std::cout << "bTilde" << std::endl;
	//std::cout << bTilde << std::endl;
			

			
		
		/**** RECREATE INSTANCE ****/
		std::cout << "recreating instance" << std::endl;
			
				
			/**** SHAPE CREATION ****/
			std::cout << "creating the shape" << std::endl;
						
				// new shape parameters
				vnl_vector<double> bShapeTilde;
				bShapeTilde = (1/w) * bTilde.extract(_nOfModes, 0);
	//std::cout << "bShapeTilde" << std::endl;
	//std::cout << bShapeTilde << std::endl;

				// allocating
				FieldType::Pointer x = FieldType::New();
				x->SetRegions(averageVFreader->GetOutput()->GetRequestedRegion());
				x->CopyInformation(averageVFreader->GetOutput());
				x->Allocate();
				x->FillBuffer(0.0);				

				// for each mode
				for (int a=0; a<_nOfModes; a++){
							
					// load the shape eigenvector
					QString modeFileName = _shapeEVectorsFileNames[0];
					if (modeFileName.lastIndexOf("/") == -1)
						modeFileName.remove(modeFileName.lastIndexOf("\\")+1,modeFileName.size()-1);
					else 
						modeFileName.remove(modeFileName.lastIndexOf("/")+1,modeFileName.size()-1);
					modeFileName.append("PC");
					modeFileName.append(QString("%1").arg(a+1));
					modeFileName.append("Shape.mhd");
					//std::cout << modeFileName.toAscii().data() << std::endl;
					FieldReaderType::Pointer shapeEVectorReader = FieldReaderType::New();
					shapeEVectorReader->SetFileName(modeFileName.toAscii().data());
					shapeEVectorReader->Update();
	//FieldRegionIteratorType p (shapeEVectorReader->GetOutput(), shapeEVectorReader->GetOutput()->GetRequestedRegion());
	//for (p.GoToBegin(); !p.IsAtEnd(); ++p)
	//	std::cout << p.Get() << std::endl;

					// multiply by the new shape parameters
					//std::cout << "multiply by bShapeTilde" << std::endl;
					FieldRegionIteratorType it6 (shapeEVectorReader->GetOutput(), shapeEVectorReader->GetOutput()->GetRequestedRegion());
					FieldRegionIteratorType it7 (x, x->GetRequestedRegion());

					for (it6.GoToBegin(), it7.GoToBegin(); !it6.IsAtEnd(); ++it6, ++it7) {
						FieldType::PixelType temp;
						temp[0] = it7.Get()[0] + it6.Get()[0] * bShapeTilde(a);
						temp[1] = it7.Get()[1] + it6.Get()[1] * bShapeTilde(a);
						temp[2] = it7.Get()[2] + it6.Get()[2] * bShapeTilde(a);
						it7.Set(temp);
	//std::cout << it7.Get()[0] << ' ' << it6.Get()[0] << ' ' << bShapeTilde(a) << std::endl;
	//std::cout << it7.Get()[1] << ' ' << it6.Get()[1] << ' ' << bShapeTilde(a) << std::endl;
	//std::cout << it7.Get()[2] << ' ' << it6.Get()[2] << ' ' << bShapeTilde(a) << std::endl;	
					}
				}
				
				// sum up the average velocity field
				std::cout << "summing up the average" << std::endl;
				FieldRegionIteratorType it8 (x, x->GetRequestedRegion());
				FieldRegionIteratorType it9 (averageVFreader->GetOutput(), averageVFreader->GetOutput()->GetRequestedRegion() );
				for (it8.GoToBegin(), it9.GoToBegin(); !it8.IsAtEnd(); ++it8, ++it9){
					it8.Set (it8.Get() + it9.Get());
				}
	//FieldRegionIteratorType h (x, x->GetRequestedRegion());
	//for (h.GoToBegin(); !h.IsAtEnd(); ++h)
	//	std::cout << h.Get() << std::endl; 

				// invert VF
				std::cout << "VF inversion" << std::endl;
				typedef itk::MultiplyByConstantImageFilter< FieldType, double, FieldType > MultiplyFilter;
				MultiplyFilter::Pointer multiplyFilter = MultiplyFilter::New();
				multiplyFilter->SetInput(x);
				multiplyFilter->SetConstant(-1.0);
				multiplyFilter->Update();

				// write inverted VF
				QString fileName = _outputFolder;
				QString temp2 = _instanceShapeFileNames[i];
				if (temp2.lastIndexOf("/") == -1){
					temp2.remove(temp2.lastIndexOf("\\")+6,temp2.size()-1);
					temp2.remove(0, temp2.lastIndexOf("\\")+1);
				}
				else {
					temp2.remove(temp2.lastIndexOf("/")+6,temp2.size()-1);
					temp2.remove(0, temp2.lastIndexOf("/")+1);
				}
				
				fileName.append("\\");
				fileName.append(temp2);
				fileName.append("_");
				fileName.append(QString("%1").arg(_nOfModes)); 
				fileName.append ("_modes_vf.mhd");
				std::cout << fileName.toAscii().data() << std::endl;
				FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
				fieldWriter->SetFileName(fileName.toAscii().data());
				fieldWriter->SetInput(multiplyFilter->GetOutput());
				fieldWriter->Update();
				
				// VF to DVF 
				std::cout << "DVF creation of recreated and original VF" << std::endl;
				typedef itk::ExponentialDeformationFieldImageFilter< FieldType,FieldType > ExponentialFieldFilterType;
				ExponentialFieldFilterType::Pointer exponentiatorRecreated = ExponentialFieldFilterType::New();
				exponentiatorRecreated->SetInput(multiplyFilter->GetOutput());
				exponentiatorRecreated->Update();

				// inverting the original VF and tranforming it to a DVF
				MultiplyFilter::Pointer multiplyFilterOriginal = MultiplyFilter::New();
				multiplyFilterOriginal->SetInput(fieldToRecreateReader->GetOutput());
				multiplyFilterOriginal->SetConstant(-1);
				multiplyFilterOriginal->Update();
				ExponentialFieldFilterType::Pointer exponentiatorOriginal= ExponentialFieldFilterType::New();
				exponentiatorOriginal->SetInput(multiplyFilterOriginal->GetOutput());
				exponentiatorOriginal->Update();
				
								

			/**** INTENSITY CREATION ****/
			std::cout << "creating the intensity" << std::endl;

						
				// combined pca members
				vnl_vector<double> bItilde;
				bItilde = bTilde.extract(_nOfModes, nOfOriginalImages);
							
				// image to create
				ImageType::Pointer g = ImageType::New();
				g->SetRegions( averageIntensityReader->GetOutput()->GetRequestedRegion() );
				g->CopyInformation( averageIntensityReader->GetOutput() );
				g->Allocate();
				g->FillBuffer(0.0);

				// for each intensity eigenvector
				for (int a=0; a<_nOfModes; a++){

					// load the intensity eigenvector
					QString modeFileName = _intensityEVectorsFileNames[0];
					if (modeFileName.lastIndexOf("/") == -1)
						modeFileName.remove(modeFileName.lastIndexOf("\\")+1,modeFileName.size()-1);
					else 
						modeFileName.remove(modeFileName.lastIndexOf("/")+1,modeFileName.size()-1);
					modeFileName.append("PC");
					modeFileName.append(QString("%1").arg(a+1));
					modeFileName.append("Intensity.mhd");
					//std::cout << modeFileName.toAscii().data() << std::endl;
					ImageReaderType::Pointer intensityEVectorReader = ImageReaderType::New();
					intensityEVectorReader->SetFileName(modeFileName.toAscii().data());
					intensityEVectorReader->Update();

					// multiply by the new intensity parameter
					ImageRegionIteratorType it10 (intensityEVectorReader->GetOutput(), intensityEVectorReader->GetOutput()->GetRequestedRegion());
					ImageRegionIteratorType it11 (g, g->GetRequestedRegion());

					for (it10.GoToBegin(), it11.GoToBegin(); !it10.IsAtEnd(); ++it10, ++it11) {
						ImageType::PixelType temp;
						temp = it11.Get() + it10.Get() * bItilde(a);
						it11.Set(temp);
					}
				}

				// sum up the average intensity
				std::cout << "sum up the average intensity" << std::endl;
				ImageRegionIteratorType it12 (g, g->GetRequestedRegion());
				ImageRegionIteratorType it13 (averageIntensityReader->GetOutput(), averageIntensityReader->GetOutput()->GetRequestedRegion() );
				for (it12.GoToBegin(), it13.GoToBegin(); !it12.IsAtEnd(); ++it12, ++it13){
					it12.Set (it12.Get() + it13.Get());
				}
				
				// write the new intensity (in the reference shape)
				ImageWriterType::Pointer imageWriter = ImageWriterType::New();
				fileName.replace(QString("vf"), QString("intensity"));
				std::cout << fileName.toAscii().data() << std::endl;
				imageWriter->SetFileName( fileName.toAscii().data());
				imageWriter->SetInput(g);
				imageWriter->Update();
							

			
			/**** COMPOSING SHAPE AND INTENSITY ****/
							
				// warping intensity to shape
				typedef itk::WarpImageFilter< ImageType, ImageType, FieldType> WarpImageFilterType;
				WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();
				typedef WarpImageFilterType::CoordRepType CoordRepType;
				typedef itk::LinearInterpolateImageFunction< ImageType,CoordRepType > InterpolatorType;
				InterpolatorType::Pointer interpolator = InterpolatorType::New();
				
				warpImageFilter->SetInput(g);
				warpImageFilter->SetDeformationField(exponentiatorRecreated->GetOutput());
				warpImageFilter->SetInterpolator(interpolator);
				warpImageFilter->SetOutputSpacing(g->GetSpacing());
				warpImageFilter->SetOutputOrigin(g->GetOrigin());
				warpImageFilter->SetEdgePaddingValue(-1024);
				warpImageFilter->Update();
			
				// writing recreated instance
				std::cout << "writing image" << std::endl;
				fileName.replace(QString("intensity"), QString("instance"));
				imageWriter->SetFileName( fileName.toAscii().data());
				imageWriter->SetInput( warpImageFilter->GetOutput() );
				imageWriter->Update();
			
				
			
			/**** STATISTICS ****/
				std::cout << "computing statistics" << std::endl;
				StatisticsDistanceCalculator* distanceCalculator = new StatisticsDistanceCalculator;
					
				// parameters - euclidean distance
				distanceCalculator->SetVectorOne(b);
				distanceCalculator->SetVectorTwo(bTilde);
				distanceCalculator->CalculateEuclideanDistance();
				vnl_vector<double> parameterDistance = distanceCalculator->GetDistanceVector();
				
				QString outputFileName = _outputFolder;
				outputFileName.append("\\");
				outputFileName.append(temp2);
				outputFileName.append("_parameterDistances_mode_");
				outputFileName.append(QString("%1").arg(_nOfModes)); 
				outputFileName.append(".txt");
				VnlWriterVector* writerVector = new VnlWriterVector;
				writerVector->SetFileName(outputFileName);
				writerVector->SetVnlVector(parameterDistance);
				writerVector->Update();

				std::cout << "parameter difference average: " << parameterDistance.mean() << std::endl;

				
				// shape - euclidean distance - shape (! on DVF to be in the image space)
				// new recreated DVF
				FieldRegionIteratorType it14 (exponentiatorRecreated->GetOutput(), exponentiatorRecreated->GetOutput()->GetRequestedRegion());
				FieldRegionIteratorType it15 (exponentiatorOriginal->GetOutput(), exponentiatorOriginal->GetOutput()->GetRequestedRegion());
	
				// calculating DVF differences
				FieldType::SizeType sizeField = fieldToRecreateReader->GetOutput()->GetBufferedRegion().GetSize();
				vnl_vector<double> shapeDistance;
				shapeDistance.set_size(sizeField[0] * sizeField[1] * sizeField[2]);
				int a=0;

				for ( it14.GoToBegin(), it15.GoToBegin(); !it14.IsAtEnd(); ++it14, ++it15){
					
					shapeDistance(a) = std::sqrt(
									(it14.Get()[0]-it15.Get()[0]) * (it14.Get()[0]-it15.Get()[0]) +
									(it14.Get()[1]-it15.Get()[1]) * (it14.Get()[1]-it15.Get()[1]) +
									(it14.Get()[2]-it15.Get()[2]) * (it14.Get()[2]-it15.Get()[2]));
					a++;
				}

				outputFileName.replace(QString("_parameterDistances_mode_"), QString("_shapeDistances_mode_"));
				writerVector->SetFileName(outputFileName);
				writerVector->SetVnlVector(shapeDistance);
				writerVector->Update();
				std::cout << "shape difference average: " << shapeDistance.mean() << std::endl;

				/*
				// compare VF before inversion
				FieldRegionIteratorType it20 (fieldToRecreateReader->GetOutput(), fieldToRecreateReader->GetOutput()->GetRequestedRegion());
				FieldRegionIteratorType it21 (x, x->GetRequestedRegion());
				vnl_vector<double> vfDistance;
				vfDistance.set_size(sizeField[0] * sizeField[1] * sizeField[2]);
				a=0;

				for ( it20.GoToBegin(), it21.GoToBegin(); !it20.IsAtEnd(); ++it20, ++it21){
					
					vfDistance(a) = std::sqrt(
									(it20.Get()[0]-it21.Get()[0]) * (it20.Get()[0]-it21.Get()[0]) +
									(it20.Get()[1]-it21.Get()[1]) * (it20.Get()[1]-it21.Get()[1]) +
									(it20.Get()[2]-it21.Get()[2]) * (it20.Get()[2]-it21.Get()[2]));
					a++;
				}	
				std::cout << "vf difference average: " << vfDistance.mean() << std::endl;
				*/

	 
				// intensity - euclidean distance 
				ImageRegionIteratorType it16 (imageToRecreateReader->GetOutput(), imageToRecreateReader->GetOutput()->GetRequestedRegion());
				ImageRegionIteratorType it17 (g, g->GetRequestedRegion());
				ImageType::SizeType sizeImage = imageToRecreateReader->GetOutput()->GetBufferedRegion().GetSize();
				vnl_vector<double> intensityDistance;
				intensityDistance.set_size(sizeImage[0] * sizeImage[1] * sizeImage[2]);

				a=0;
				for ( it16.GoToBegin(), it17.GoToBegin(); !it16.IsAtEnd(); ++it16, ++it17){
					
					intensityDistance(a) = std::sqrt((it16.Get()-it17.Get()) * (it16.Get()-it17.Get()));
					a++;
				}

				outputFileName.replace(QString("_shapeDistances_mode_"), QString("_intensityDistances_mode_"));
				writerVector->SetFileName(outputFileName);
				writerVector->SetVnlVector(intensityDistance);
				writerVector->Update();
				
				std::cout << "intensity difference average: " << intensityDistance.mean() << std::endl;

				averageParameterDistances(i,m) = parameterDistance.mean();
				averageShapeDistances(i,m) = shapeDistance.mean();
				averageIntensityDistances(i,m) = intensityDistance.mean() ;

	//std::cout << outputFileName.toAscii().data() << std::endl;

	//std::cout << "original image" << std::endl; 
	//for ( it16.GoToBegin(); !it16.IsAtEnd(); ++it16)
	//	std::cout << it16.Get() << std::endl;
	//std::cout << "recreated image" << std::endl; 
	//for ( it17.GoToBegin(); !it17.IsAtEnd(); ++it17)
	//	std::cout << it17.Get() << std::endl;
	//std::cout << "intensityDistance" << std::endl;
	//std::cout << intensityDistance << std::endl;

				
				double end = clock();
				double total = (end-start)/CLOCKS_PER_SEC;
				std::cout << "Computation time: " << total << " sec. (" << total/60 << " min.)" << std::endl;

						
		}

	}


		/**** GLOBAL STATISTICS ****/ 
	
		VnlWriterMatrix* writer = new VnlWriterMatrix;

		writer->SetVnlMatrix(averageParameterDistances);
		_outputFolder.append("\\averageParameterDistances.txt");
		writer->SetFileName(_outputFolder);
		writer->MatrixShapeUpdate();

		writer->SetVnlMatrix(averageShapeDistances);
		_outputFolder.replace(QString("averageParameterDistances"), QString ("averageShapeDistances"));
		writer->SetFileName(_outputFolder);
		writer->MatrixShapeUpdate();

		writer->SetVnlMatrix(averageIntensityDistances);
		_outputFolder.replace(QString("averageShapeDistances"), QString ("averageIntensityDistances"));
		writer->SetFileName(_outputFolder);
		writer->MatrixShapeUpdate();

		delete writer;

		std::cout << "global averages written" << std::endl;

	}
	
	// private - test
	void PCAImages::RecreateShapePca(){

		// load the velocity field to recreate
		std::cout << "loading the shape to recreate" << std::endl;
		FieldReaderType::Pointer fieldToRecreateReader = FieldReaderType::New();
		std::cout << _instanceShapeFileNames[0].toAscii().data() << std::endl;
		fieldToRecreateReader->SetFileName(_instanceShapeFileNames[0].toAscii().data());
		fieldToRecreateReader->Update();
//FieldRegionIteratorType q (fieldToRecreateReader->GetOutput(), fieldToRecreateReader->GetOutput()->GetRequestedRegion());
//for (q.GoToBegin(); !q.IsAtEnd(); ++q)
//	std::cout << q.Get() << std::endl;

		
		// load average velocity field
		std::cout << "loading the average shape" << std::endl;
		FieldReaderType::Pointer averageVFreader = FieldReaderType::New();
		std::cout << _shapeAverageFileName.toAscii().data() << std::endl;
		averageVFreader->SetFileName(_shapeAverageFileName.toAscii().data());
		averageVFreader->Update();
//FieldRegionIteratorType w (averageVFreader->GetOutput(), averageVFreader->GetOutput()->GetRequestedRegion());
//for (w.GoToBegin(); !w.IsAtEnd(); ++w)
//	std::cout << w.Get() << std::endl;


		// subtract meanVelField from instance
std::cout << "subtracting the average" << std::endl;
		typedef itk::SubtractImageFilter< FieldType, FieldType, FieldType> SubtractFieldFilterType;
		SubtractFieldFilterType::Pointer subtractFieldFilter = SubtractFieldFilterType::New();
		subtractFieldFilter->SetInput1(fieldToRecreateReader->GetOutput());
		subtractFieldFilter->SetInput2(averageVFreader->GetOutput());
		subtractFieldFilter->Update();
//FieldRegionIteratorType e (subtractFieldFilter->GetOutput(), subtractFieldFilter->GetOutput()->GetRequestedRegion());
//for (e.GoToBegin(); !e.IsAtEnd(); ++e)
//	std::cout << e.Get() << std::endl;

		
		vnl_vector<double> bs;
		bs.set_size(6);

		// calculate the parameters
		for(int j = 0; j < 6; ++j) {

			// empty image (sum of the weighted combined modes)
			FieldType::Pointer fiBsProduct = FieldType::New();
			fiBsProduct->SetRegions( averageVFreader->GetOutput()->GetRequestedRegion() );
			fiBsProduct->CopyInformation( averageVFreader->GetOutput() );
			fiBsProduct->Allocate();
			fiBsProduct->FillBuffer(0.0);
			
			// read the i-th shape eigenvector
			FieldReaderType::Pointer shapeEVectorReader = FieldReaderType::New();
			QString temp = ("PC");
			temp.append(QString("%1").arg(j+1));
			temp.append("Shape.mhd");
			for (int y=0; y<_shapeEVectorsFileNames.size(); y++){
				if (_shapeEVectorsFileNames[y].contains(temp)){
					shapeEVectorReader->SetFileName(_shapeEVectorsFileNames[y].toAscii().data());
					std::cout << _shapeEVectorsFileNames[y].toAscii().data() << std::endl;
				}
			}
			shapeEVectorReader->Update();
//FieldRegionIteratorType r (shapeEVectorReader->GetOutput(), shapeEVectorReader->GetOutput()->GetRequestedRegion());
//for (r.GoToBegin(); !r.IsAtEnd(); ++r)
//	std::cout << r.Get() << std::endl;

			// multiply by the difference
			FieldRegionIteratorType it0 (shapeEVectorReader->GetOutput(), shapeEVectorReader->GetOutput()->GetRequestedRegion());
			FieldRegionIteratorType it1 (subtractFieldFilter->GetOutput(), subtractFieldFilter->GetOutput()->GetRequestedRegion());
			FieldRegionIteratorType it2 (fiBsProduct, fiBsProduct->GetRequestedRegion());
			
			for (it0.GoToBegin(), it1.GoToBegin(), it2.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1, ++it2) {
				FieldType::PixelType value;
				value[0] = it0.Get()[0] * it1.Get()[0] + it2.Get()[0];
				value[1] = it0.Get()[1] * it1.Get()[1] + it2.Get()[1];
				value[2] = it0.Get()[2] * it1.Get()[2] + it2.Get()[2];
				it2.Set(value);
			}

			// sum of the three parameters
			double sum0 = 0.0;
			for (it2.GoToBegin(); !it2.IsAtEnd(); ++it2) {
				sum0 += it2.Get()[0];
				sum0 += it2.Get()[1];
				sum0 += it2.Get()[2];
			}
						
			bs(j) = sum0; //bs(j+nOfOriginalImages) = sum0; bs(j+2*nOfOriginalImages) = sum0;
		}

		std::cout << "bs" << std::endl;
		std::cout << bs << std::endl;



		// recreate the shape
		// allocating
		FieldType::Pointer x = FieldType::New();
		x->SetRegions(averageVFreader->GetOutput()->GetRequestedRegion());
		x->CopyInformation(averageVFreader->GetOutput());
		x->Allocate();
		x->FillBuffer(0.0);				

		// for each mode
		for (int a=0; a<_nOfModes; a++){
					
			// load the shape eigenvectors
			QString modeFileName = _shapeEVectorsFileNames[0];
			if (modeFileName.lastIndexOf("/") == -1)
				modeFileName.remove(modeFileName.lastIndexOf("\\")+1,modeFileName.size()-1);
			else 
				modeFileName.remove(modeFileName.lastIndexOf("/")+1,modeFileName.size()-1);
			modeFileName.append("PC");
			modeFileName.append(QString("%1").arg(a+1));
			modeFileName.append("Shape.mhd");
			std::cout << modeFileName.toAscii().data() << std::endl;
			FieldReaderType::Pointer shapeEVectorReader = FieldReaderType::New();
			shapeEVectorReader->SetFileName(modeFileName.toAscii().data());
			shapeEVectorReader->Update();

			// multiply by the new shape parameters
			FieldRegionIteratorType it6 (shapeEVectorReader->GetOutput(), shapeEVectorReader->GetOutput()->GetRequestedRegion());
			FieldRegionIteratorType it7 (x, x->GetRequestedRegion());

			for (it6.GoToBegin(), it7.GoToBegin(); !it6.IsAtEnd(); ++it6, ++it7) {
				FieldType::PixelType temp;
				temp[0] = it7.Get()[0] + it6.Get()[0] * bs(a);
				temp[1] = it7.Get()[1] + it6.Get()[1] * bs(a);
				temp[2] = it7.Get()[2] + it6.Get()[2] * bs(a);
				it7.Set(temp);
//std::cout << it7.Get()[0] << ' ' << it6.Get()[0] << ' ' << bs(a) << std::endl;
//std::cout << it7.Get()[1] << ' ' << it6.Get()[1] << ' ' << bs(a) << std::endl;
//std::cout << it7.Get()[2] << ' ' << it6.Get()[2] << ' ' << bs(a) << std::endl;	
			}
		}
		
		// sum up the average velocity field
		std::cout << "summing up the average" << std::endl;
		FieldRegionIteratorType it8 (x, x->GetRequestedRegion());
		FieldRegionIteratorType it9 (averageVFreader->GetOutput(), averageVFreader->GetOutput()->GetRequestedRegion() );
		for (it8.GoToBegin(), it9.GoToBegin(); !it8.IsAtEnd(); ++it8, ++it9){
			it8.Set (it8.Get() + it9.Get());
		}

//FieldRegionIteratorType t (x, x->GetRequestedRegion());
//for (t.GoToBegin(); !t.IsAtEnd(); ++t)
//	std::cout << t.Get() << std::endl;

		// write vf
		QString fileName = _outputFolder;
		
		if (_instanceShapeFileNames[0].lastIndexOf("/") == -1){
			_instanceShapeFileNames[0].remove(_instanceShapeFileNames[0].lastIndexOf("\\")+6,_instanceShapeFileNames[0].size()-1);
			_instanceShapeFileNames[0].remove(0, _instanceShapeFileNames[0].lastIndexOf("\\")+1);
		}
		else {
			_instanceShapeFileNames[0].remove(_instanceShapeFileNames[0].lastIndexOf("/")+6,_instanceShapeFileNames[0].size()-1);
			_instanceShapeFileNames[0].remove(0, _instanceShapeFileNames[0].lastIndexOf("/")+1);
		}
		
		fileName.append("\\");
		fileName.append(_instanceShapeFileNames[0]);
		fileName.append("_recreated_vf_");
		fileName.append(QString("%1").arg(_nOfModes)); 
		fileName.append (".mhd");
		std::cout << fileName.toAscii().data() << std::endl;
		FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
		fieldWriter->SetFileName(fileName.toAscii().data());
		fieldWriter->SetInput(x);
		//fieldWriter->Update();

		FieldRegionIteratorType it20 (fieldToRecreateReader->GetOutput(), fieldToRecreateReader->GetOutput()->GetRequestedRegion());
		FieldRegionIteratorType it21 (x, x->GetRequestedRegion());
		vnl_vector<double> vfDistance;
		FieldType::SizeType sizeField = fieldToRecreateReader->GetOutput()->GetBufferedRegion().GetSize();
		vfDistance.set_size(sizeField[0] * sizeField[1] * sizeField[2]);
		int a=0;
		for ( it20.GoToBegin(), it21.GoToBegin(); !it20.IsAtEnd(); ++it20, ++it21){
			
			vfDistance(a) = std::sqrt(
							(it20.Get()[0]-it21.Get()[0]) * (it20.Get()[0]-it21.Get()[0]) +
							(it20.Get()[1]-it21.Get()[1]) * (it20.Get()[1]-it21.Get()[1]) +
							(it20.Get()[2]-it21.Get()[2]) * (it20.Get()[2]-it21.Get()[2]));
			a++;
		}	

		std::cout << "vf difference average: " << vfDistance.mean() << std::endl;
		



	}
}

/*
			// extract parameters (the same as _nOfModes)
			vnl_vector<double> parameterS;
			parameterS.set_size(_nOfModes * 3);
			for (int j=0; j<_nOfModes * 3; j++){
				parameterS(j) = bs(j);
			}
			for (int j=0; j<_nOfModes; j++){
				parameterS(j) = bs(j);
				std::cout << "eigenvalue " << j+1 << " weight: " << parameterS(j)/std::sqrt(combinedEValues(j)) << std::endl;
			}
			std::cout << "shape parameters: " <<  parameterS << std::endl;

			// empty image (sum of the weighted combined modes)
			FieldType::Pointer x = FieldType::New();
			x->SetRegions( averageVFreader->GetOutput()->GetRequestedRegion() );
			x->CopyInformation( averageVFreader->GetOutput() );
			x->Allocate();
			x->FillBuffer(0.0);

			// for each mode
			for (int n=0; n<_nOfModes; n++){
				std::cout << "mode: " << n+1 << std::endl;
				
				// for each shape eigenvector
				for (int a=0; a<nOfOriginalImages; a++){
							
					// load the a-th shape eigenvector
					QString modeFileName = _shapeEVectorsFileNames[0];
					if (modeFileName.lastIndexOf("/") == -1)
						modeFileName.remove(modeFileName.lastIndexOf("\\")+1,modeFileName.size()-1);
					else 
						modeFileName.remove(modeFileName.lastIndexOf("/")+1,modeFileName.size()-1);
					modeFileName.append("PC");
					modeFileName.append(QString("%1").arg(a+1));
					modeFileName.append("Shape.mhd");
					std::cout << modeFileName.toAscii().data() << std::endl;
					
					ImageHandler* shapeEVectorReader = new ImageHandler;
					shapeEVectorReader->SetFieldFileName(modeFileName);
					shapeEVectorReader->FieldReaderUpdate();
					
					// multiply the a-th shape eigenvector by 1/w
					std::cout << "multiply by 1/w: " << 1.0/w << std::endl;
					typedef itk::MultiplyByConstantImageFilter< FieldType, double, FieldType > MultiplyFilter;
					MultiplyFilter::Pointer multiplyFilter = MultiplyFilter::New();
					multiplyFilter->SetInput(shapeEVectorReader->GetField());
					multiplyFilter->SetConstant(1.0/w);
					multiplyFilter->Update();

					typedef itk::ImageRegionIterator< FieldType > RegionIteratorType;
					RegionIteratorType it10 (multiplyFilter->GetOutput(), multiplyFilter->GetOutput()->GetRequestedRegion());
					
					// multiply the result by the a-th voxel position of the n-th combined eigenvector
					std::cout << "multiply by the combined eigenvector values" << std::endl;
					
					// the n-th combined eigenvector 
					float firstElementCombinedEV = combinedEigenVectors(a, n);
					float secondElementCombinedEV = combinedEigenVectors(a + nOfOriginalImages, n);
					float thirdElementCombinedEV = combinedEigenVectors(a + 2*nOfOriginalImages, n);
					//std::cout << n+1 << " combined parameter: " << firstElementCombinedEV << ' ' << secondElementCombinedEV << ' ' << thirdElementCombinedEV << std::endl;
					
					// multiply the a-th eigenvector by the a-th position of the n-th combined eigenvector
					FieldRegionIteratorType it0 (multiplyFilter->GetOutput(), multiplyFilter->GetOutput()->GetRequestedRegion());
					FieldRegionIteratorType it1 (x, x->GetRequestedRegion());
	
					for (it0.GoToBegin(), it1.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1) {
												
						// multiply	(3 values since they are vector images)			
						float value0 = it0.Get()[0] * firstElementCombinedEV;
						float value1 = it0.Get()[1] * secondElementCombinedEV;
						float value2 = it0.Get()[2] * thirdElementCombinedEV;
	
						// multiply by the weighted eigenvalue
						value0 *= parameterS(n);
						value1 *= parameterS(n);
						value2 *= parameterS(n);
								
						// update the sum of velocity fields
						FieldType::PixelType temp;
						temp[0] = it1.Get()[0] + value0; temp[1] = it1.Get()[1] + value1; temp[2] = it1.Get()[2] + value2;
						it1.Set(temp);
					}
					// cleaning
					delete shapeEVectorReader;
				}
			}
			
			
			
			FieldRegionIteratorType it20 (x, x->GetRequestedRegion());
			
			std::cout << "without optimization" << std::endl;
			for ( it20.GoToBegin(); !it20.IsAtEnd(); ++it20){
				std::cout << it20.Get() << std::endl;
			}
			*/


/*

			// extract parameters (bg contains nOfOriginalImages parameters, parameterG the needed _nOfModes)
			vnl_vector<double> parameterG;
			parameterG.set_size(_nOfModes);
			for (int j=0; j<_nOfModes; j++){
				parameterG(j) = bg(j);
				std::cout << "eigenvalue " << j+1 << "parameterG: " << parameterG << " weight: " << parameterG(j)/std::sqrt(combinedEValues(j)) << std::endl;
			}

			// empty image (sum of the weighted combined modes)
			ImageType::Pointer g = ImageType::New();
			g->SetRegions( averageIntensityReader->GetOutput()->GetRequestedRegion() );
			g->CopyInformation( averageIntensityReader->GetOutput() );
			g->Allocate();
			g->FillBuffer(0.0);

			// for each mode
			for (int n=0; n<_nOfModes; n++){
				std::cout << "mode: " << n+1 << std::endl;

				// for each intensity eigenvector
				for (int a=0; a<nOfOriginalImages; a++){

					// load the a-th intensity eigenvector
					QString modeFileName = _intensityEVectorsFileNames[0];
					if (modeFileName.lastIndexOf("/") == -1)
						modeFileName.remove(modeFileName.lastIndexOf("\\")+1,modeFileName.size()-1);
					else 
						modeFileName.remove(modeFileName.lastIndexOf("/")+1,modeFileName.size()-1);
					modeFileName.append("PC");
					modeFileName.append(QString("%1").arg(a+1));
					modeFileName.append("Intensity.mhd");
					std::cout << modeFileName.toAscii().data() << std::endl;
					
					ImageReaderType::Pointer intensityEVectorReader = ImageReaderType::New();
					intensityEVectorReader->SetFileName(modeFileName.toAscii().data());
					intensityEVectorReader->Update();

					// multiply the result by the a-th voxel position of the n-th combined eigenvector
					std::cout << "multiply by the combined eigenvector values" << std::endl;
					
					// the n-th combined eigenvector 
					float fourthElementCombinedEV = combinedEigenVectors(a + 3*nOfOriginalImages, n);
					std::cout << n+1 << " combined parameter: " << fourthElementCombinedEV << std::endl;
					std::cout << "parameter c: " << parameterG(n) << std::endl;
					
					ImageRegionIteratorType it0 (intensityEVectorReader->GetOutput(), intensityEVectorReader->GetOutput()->GetRequestedRegion());
					ImageRegionIteratorType it1 (g, g->GetRequestedRegion());
	
					for (it0.GoToBegin(), it1.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1) {
												
						// multiply	(3 values since they are vector images)			
						float value4 = it0.Get() * fourthElementCombinedEV;
						//std::cout << value4 << ' ' << it0.Get() << ' ' << fourthElementCombinedEV << std::endl;

						// multiply by the weighted eigenvalue
						value4 *= parameterG(n);
						//std::cout << value4 << ' ' << parameterG(n) <<std::endl;

						// update the sum of images
						ImageType::PixelType temp;
						temp = it1.Get() + value4;
						///std::cout << temp << ' ' << it1.Get() << ' ' << value4 << std::endl;
						it1.Set(temp);
						//std::cout << it1.Get() << std::endl; 
					}
				}
			}
			
			std::cout << "without optimization" << std::endl;
			typedef itk::ImageRegionIterator< ImageType > IteratorType;
			IteratorType it9 (g, g->GetRequestedRegion());
			for ( it9.GoToBegin(); !it9.IsAtEnd(); ++it9){
				std::cout << it9.Get() << std::endl;
			}
*/

