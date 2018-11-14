/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PCAImageGeneralization.h>

#include <itkImageFileReader.h>
#include <itkSubtractImageFilter.h>
#include <itkImageRegionIterator.h>

namespace pca{

	// constructor
	PCAImageGeneralization::PCAImageGeneralization(){

		_distanceVector = vtkDoubleArray::New();
		_distanceAverage = 0;
		_distanceStandardDeviation = 0;
		_distanceStandardError = 0;

	}

	// destructor
	PCAImageGeneralization::~PCAImageGeneralization(){
	}
	
	// overwritten virtual function
	void PCAImageGeneralization::Update(){

		for (int i=0; i<_datasetFileNames.size(); i++){
			
			// load the image dataset 
			typedef itk::ImageFileReader< ImageType > ReaderType;
			ReaderType::Pointer datasetReader = ReaderType::New();
			datasetReader->SetFileName(_datasetFileNames[i].ascii());
			//std::cout << "dataset: " << _datasetFileNames[i].ascii() << std::endl;
			datasetReader->Update();

			double minValue;
			
			for (int y=0; y<_instanceFileNames.size(); y++){
				
				// load the instance
				ReaderType::Pointer instanceReader = ReaderType::New();
				instanceReader->SetFileName(_instanceFileNames[y].ascii());
				//std::cout << "instance: " << _instanceFileNames[y].ascii() << std::endl;
				instanceReader->Update();

				// subtract the two images
				typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType > SubtractImageFilterType;
				SubtractImageFilterType::Pointer subtractor = SubtractImageFilterType::New();
				subtractor->SetInput1(datasetReader->GetOutput());
				subtractor->SetInput2(instanceReader->GetOutput());
				subtractor->Update();

				// sum up the differences 
				typedef itk::ImageRegionIterator< ImageType > ImageRegionIteratorType;
				ImageRegionIteratorType it (subtractor->GetOutput(), subtractor->GetOutput()->GetRequestedRegion());
				
				float sum = 0;
				for (it.GoToBegin(); !it.IsAtEnd(); ++it){
					sum += abs(it.Get());
				}
				//std::cout << "sum: " << sum << std::endl;

				// looking for the minimum distance (closest bone)
				if (y==0)
					minValue = sum;
				if (sum<minValue)
					minValue = sum;
				
			}

//// divide by the number of voxels!!!!! ///////////////		

			_distanceVector->InsertNextValue(minValue);
		
		
		// average
		_distanceAverage = _distanceAverage + (1.0/_datasetFileNames.size()) * minValue;
		//std::cout << "distanceAverage: " << _distanceAverage << std::endl;
		}

	// standard deviation
	for (int a=0; a<_datasetFileNames.size(); a++){
		_distanceStandardDeviation = _distanceStandardDeviation + (_distanceVector->GetValue(a) - _distanceAverage)*(_distanceVector->GetValue(a) - _distanceAverage);
	}
	_distanceStandardDeviation = _distanceStandardDeviation * (1.0/(_datasetFileNames.size()-1));
	_distanceStandardDeviation = std::sqrt (_distanceStandardDeviation);
	//std::cout << "_distanceStandardDeviation: " << _distanceStandardDeviation << std::endl;

	// standard error
	_distanceStandardError = _distanceStandardDeviation / std::sqrt(double(_datasetFileNames.size()));
	//std::cout << "_distanceStandardError: " << _distanceStandardError << std::endl;
	}
}
