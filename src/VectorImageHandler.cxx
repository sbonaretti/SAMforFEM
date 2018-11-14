/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <VectorImageHandler.h>

#include <itkExponentialDeformationFieldImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkVTKImageToImageFilter.h>

#include <vtkDoubleArray.h>

namespace image{

	// constructor
	VectorImageHandler::VectorImageHandler(){

	}

	// destructor
	VectorImageHandler::~VectorImageHandler(){

		_field = NULL;
		//_vtkField ->Delete();
	}

	// member functions
	
	/**************************** FIELD READER AND WRITER ****************************/
	void VectorImageHandler::MetafileReaderUpdate(){
				
		typedef itk::ImageFileReader< VectorImageHandler::FieldType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(_fieldFileName.toAscii().data());
		//std::cout << _fieldFileName.toAscii().data() << std::endl;
		try{
			reader->Update();
		}
		catch( itk::ExceptionObject & excp ){
			std::cerr << "!!! Problem reading the input file" << std::endl;
		}
		
		// field
		_field = reader->GetOutput();
		
		// field data
		const VectorImageHandler::FieldType::SpacingType& spacing = _field->GetSpacing();
		const VectorImageHandler::FieldType::PointType& origin = _field->GetOrigin();
		const VectorImageHandler::FieldType::SizeType& size = _field->GetBufferedRegion().GetSize();
		std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
		std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
		std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
		typedef itk::MinimumMaximumImageCalculator< VectorImageHandler::FieldType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_field);
		std::cout << "min intensity: " << calculator->GetMinimum() << "; max intensity: " << calculator->GetMaximum() << std::endl; 

	}


	void VectorImageHandler::MetafileWriterUpdate(){
		
		typedef itk::ImageFileWriter< VectorImageHandler::FieldType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( _fieldFileName.toAscii().data() );
		writer->SetInput( _field );
		writer->Update();
		std::cout << "metafile written" << std::endl;

	
	}	
	
	/****************************** ITK VTK CONVERSION *******************************/
	void VectorImageHandler::ITKtoVTK(){

		typedef itk::ImageToVTKImageFilter< FieldType > FieldToVTKImageFilterType;
		FieldToVTKImageFilterType::Pointer fieldToVTKImageFilter = FieldToVTKImageFilterType::New();
		fieldToVTKImageFilter->SetInput(_field);
		fieldToVTKImageFilter->Update();

		vtkImageData* vtkFieldImage = vtkImageData::New();
		vtkFieldImage->DeepCopy(fieldToVTKImageFilter->GetOutput());
		vtkFieldImage->Update();
	}

	void VectorImageHandler::ITKtoVTKtoITK(){
	
		// ITK to VTK
		typedef itk::ImageToVTKImageFilter< FieldType > FieldToVTKImageFilterType;
		FieldToVTKImageFilterType::Pointer fieldToVTKImageFilter = FieldToVTKImageFilterType::New();
		fieldToVTKImageFilter->SetInput(_field);
		fieldToVTKImageFilter->Update();

		vtkImageData* vtkFieldImage = vtkImageData::New();
		vtkFieldImage->DeepCopy(fieldToVTKImageFilter->GetOutput());
		vtkFieldImage->Update();
		
		// VTK to ITK
		typedef itk::VTKImageToImageFilter< FieldType > FieldConnectorType;
		FieldConnectorType::Pointer fieldConnector = FieldConnectorType::New();
		fieldConnector->SetInput(vtkFieldImage);
		fieldConnector->Update();
		
		_field = fieldConnector->GetOutput();
	
	}

	
	/******************************* VF DVF CONVERSION *******************************/
	void VectorImageHandler::VFtoDVFinverted(){

		// invert VF
		typedef itk::MultiplyByConstantImageFilter< FieldType, double, FieldType > MultiplyFilter;
		MultiplyFilter::Pointer multiplyFilter = MultiplyFilter::New();
		multiplyFilter->SetInput(_field);
		multiplyFilter->SetConstant(-1);
		multiplyFilter->Update();

		// VF to DVF 
		typedef itk::ExponentialDeformationFieldImageFilter< FieldType,FieldType > ExponentialFieldFilterType;
		ExponentialFieldFilterType::Pointer exponentiator = ExponentialFieldFilterType::New();
		exponentiator->SetInput(multiplyFilter->GetOutput());
		exponentiator->Update();
		_field = exponentiator->GetOutput();
	
	}

	void VectorImageHandler::VFtoDVFnotInverted(){

		// VF to DVF 
		typedef itk::ExponentialDeformationFieldImageFilter< FieldType,FieldType > ExponentialFieldFilterType;
		ExponentialFieldFilterType::Pointer exponentiator = ExponentialFieldFilterType::New();
		exponentiator->SetInput(_field);
		exponentiator->Update();
		_field = exponentiator->GetOutput();
	
	}

	void VectorImageHandler::SelectReferenceAsMinimumDistanceToAverage(){
	
		
		// create empty average image
		typedef itk::ImageFileReader< VectorImageHandler::FieldType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(_fieldFileNames[0]);
		try{
				reader->Update();
		}
		catch( itk::ExceptionObject & excp ){
			std::cerr << "!!! Problem reading the input file" << std::endl;
		}
				
		FieldType::Pointer averageImage = FieldType::New();
		averageImage->SetRegions(reader->GetOutput()->GetRequestedRegion());
		averageImage->SetSpacing(reader->GetOutput()->GetSpacing());
		averageImage->Allocate();
		FieldType::PixelType value;
		value.Fill(0.0);
		averageImage->FillBuffer(value);
		
		
		// calculate average
		std::cout << "calculate average velocity field" << std::endl;
		for (int i=0; i<_fieldFileNames.size(); i++){

			// load vf
			typedef itk::ImageFileReader< VectorImageHandler::FieldType > ReaderType;
			ReaderType::Pointer fieldReader = ReaderType::New();
			fieldReader->SetFileName(_fieldFileNames[i]);
			std::cout << _fieldFileNames[i].toAscii().data() << std::endl;
			try{
				fieldReader->Update();
			}
			catch( itk::ExceptionObject & excp ){
				std::cerr << "!!! Problem reading the input file" << std::endl;
			}

			// divide its values by the number of images
			typedef itk::MultiplyByConstantImageFilter< FieldType, double, FieldType > MultiplyFilter;
			MultiplyFilter::Pointer multiplyFilter = MultiplyFilter::New();
			multiplyFilter->SetInput(fieldReader->GetOutput());
			multiplyFilter->SetConstant(1.0/_fieldFileNames.size());
			multiplyFilter->Update();

			// sum up to the average field
			typedef itk::ImageRegionIterator< FieldType > IteratorType;
			IteratorType it1( averageImage, averageImage->GetRequestedRegion() );
			IteratorType it2( multiplyFilter->GetOutput(), multiplyFilter->GetOutput()->GetRequestedRegion() );
			for ( it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2 ){
				it1.Set( it1.Get() + it2.Get() );
			}
		}
		std::cout << std::endl;
		
		
		// calculate norms
		std::cout << "calculate distance norms" << std::endl;
		vtkDoubleArray* norms = vtkDoubleArray::New();

		for (int i=0; i<_fieldFileNames.size(); i++){

			// re-load vf
			typedef itk::ImageFileReader< VectorImageHandler::FieldType > ReaderType;
			ReaderType::Pointer fieldReader = ReaderType::New();
			fieldReader->SetFileName(_fieldFileNames[i]);
			std::cout << _fieldFileNames[i].toAscii().data() << std::endl;
			try{
				fieldReader->Update();
			}
			catch( itk::ExceptionObject & excp ){
				std::cerr << "!!! Problem reading the input file" << std::endl;
			}

			// calculate norm
			double norm = 0.0;
			typedef itk::ImageRegionIterator< FieldType > IteratorType;
			IteratorType it1( averageImage, averageImage->GetRequestedRegion() );
			IteratorType it2( fieldReader->GetOutput(), fieldReader->GetOutput()->GetRequestedRegion() );
			for ( it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2){
				norm += std::sqrt(( it1.Get() - it2.Get() ) * ( it1.Get() - it2.Get() ));
			}
			norms->InsertNextValue(norm);
			std::cout << "norm: " << norm << std::endl;
		}
		std::cout << std::endl;

		// select the minimum norm and print it out
		double range[2];
		norms->GetRange(range);
		std::cout << "min: " << range[0] << std::endl;
		for (int i=0; i<norms->GetSize(); i++){
			if (norms->GetValue(i) == range[0]){
				std::cout << "new reference: " <<  _fieldFileNames[i].toAscii().data() << std::endl;
				std::cout << "min norm: " << norms->GetValue(i);
			}
		}

	}
}