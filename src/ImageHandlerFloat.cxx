/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <ImageHandlerFloat.h>

#include <itkAffineTransform.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageToVTKImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkWarpImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkVTKImageToImageFilter.h>

#include <vtkMath.h>



namespace image {

	// constructor
	ImageHandlerFloat::ImageHandlerFloat(){

		_threshold = -1000;
		_vtkImage = vtkImageData::New();
		_edgePaddingValue = 0.0;
		_min = -90;
					
	}

	// destructor
	ImageHandlerFloat::~ImageHandlerFloat(){
		
		_image = NULL;

	}

	
	// member functions
	void ImageHandlerFloat::MetafileReaderUpdate(){
				
		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(_imageFileName.toAscii().data());
		//std::cout << _imageFileName.toAscii().data() << std::endl;
		try{
			reader->Update();
		}
		catch( itk::ExceptionObject & excp ){
			std::cerr << "!!! Problem reading the input file" << std::endl;
		}
		
		// image
		_image = reader->GetOutput();
		
		// image data
		const ImageType::SpacingType& spacing = _image->GetSpacing();
		const ImageType::PointType& origin = _image->GetOrigin();
		const ImageType::SizeType& size = _image->GetBufferedRegion().GetSize();
		std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
		std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
		std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
		typedef itk::MinimumMaximumImageCalculator< ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_image);
		calculator->Compute();
		std::cout << "min intensity: " << calculator->GetMinimum() << "; max intensity: " << calculator->GetMaximum() << std::endl; 
	
	}

	void ImageHandlerFloat::MetafileWriterUpdate(){

		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( _imageFileName.toAscii().data() );
		writer->SetInput( _image );
		writer->Update();
		std::cout << "metafile written" << std::endl;
	
	}

	void ImageHandlerFloat::PseudoCalibration(){	

		typedef itk::MinimumMaximumImageCalculator< ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_image);
		calculator->Compute();
		/*
		double zeroDensity = 0.0; double maxDensity = 1.08; // y-axis // this is rhoAsh
		double zeroHU = 0.0; double maxHU = calculator->GetMaximum(); // x-axis
		//double zeroHU = -90.0; double maxHU = calculator->GetMaximum(); // x-axis
		//double zeroHU = _min; double maxHU = calculator->GetMaximum(); // x-axis
		*/
		double zeroDensity = 0.10; double maxDensity = 1.08; // y-axis // this is rhoAsh
		double zeroHU = -90; double maxHU = 1650;

		double slope = (maxDensity - zeroDensity)/(maxHU - zeroHU); // slope
		double intercept = zeroDensity - ((zeroDensity-maxDensity)/(zeroHU-maxHU))*zeroHU;
		
		std::cout << "slope: " << slope << "  intercept: " << intercept << std::endl;
		
		typedef itk::ImageRegionIterator< ImageType> IteratorType;
		IteratorType it (_image, _image->GetRequestedRegion());

		ImageType::IndexType pixelIndex;
		pixelIndex[0] = 2; pixelIndex[1] = 2; pixelIndex[2] = 2;
		ImageType::PixelType backgroundValue = _image->GetPixel( pixelIndex );


		for (it.GoToBegin(); !it.IsAtEnd(); ++it){

			double greyLevel = it.Get();
			
			if (greyLevel > -800){
			//if (greyLevel != backgroundValue){
				greyLevel *= slope;
				greyLevel += intercept;
				greyLevel *= 1000; // multiplication by 1000 in order to keep short images (could be considered as g/mm3)
				it.Set(greyLevel);}
			else
				it.Set(-1024);
			
		}

	}	
	void ImageHandlerFloat::FindMinMarrow(){
	
		// warping the mask to the instance
		typedef itk::WarpImageFilter< ImageType,ImageType, FieldType> WarpImageFilterType;
		WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();
		typedef WarpImageFilterType::CoordRepType CoordRepType;
		typedef itk::NearestNeighborInterpolateImageFunction< ImageType, CoordRepType> InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		
		warpImageFilter->SetInput(_mask);
		warpImageFilter->SetDeformationField(_svf);
		warpImageFilter->SetInterpolator(interpolator);
		warpImageFilter->SetOutputSpacing(_mask->GetSpacing());
		warpImageFilter->SetOutputOrigin(_mask->GetOrigin());
		warpImageFilter->SetEdgePaddingValue(_edgePaddingValue);
		warpImageFilter->Update();
		_mask = warpImageFilter->GetOutput();

		// allocate a third image
		ImageType::Pointer marrowImage = ImageType::New();
		marrowImage->SetRegions( _image->GetRequestedRegion() );
		marrowImage->CopyInformation( _image );
		marrowImage->Allocate();
		marrowImage->FillBuffer(100);
		
		// region iterator
		typedef itk::ImageRegionIterator< ImageType> IteratorType;
		IteratorType it0 (_mask, _mask->GetRequestedRegion());
		IteratorType it1 (_image, _image->GetRequestedRegion());
		IteratorType it2 (marrowImage, marrowImage->GetRequestedRegion());

		/*
		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName("_mask.mhd");
		writer->SetInput( _mask );
		writer->Update();
		std::cout << "_mask.mhd written" << std::endl;

		writer->SetFileName("_image.mhd");
		writer->SetInput( _image );
		writer->Update();
		std::cout << "_image.mhd written" << std::endl;
		
		writer->SetFileName("_marrowImage.mhd");
		writer->SetInput( marrowImage );
		writer->Update();
		std::cout << "_marrowImage.mhd written" << std::endl;
		*/

		// extract the marrow
		for (it0.GoToBegin(), it1.GoToBegin(), it2.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1, ++it2){
			if ((it0.Get() < 110.0) && (it0.Get() > 90.0)){
				it2.Set(it1.Get());
			}
		}
		
		// extract the minimum
		typedef itk::MinimumMaximumImageCalculator< ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(marrowImage);
		calculator->Compute();
		_min = calculator->GetMinimum();
		std::cout << "min intensity: " << _min << "; max intensity: " << calculator->GetMaximum() << std::endl; 
	
		/*
		writer->SetFileName("marrow.mhd");
		writer->SetInput( marrowImage );
		writer->Update();
		std::cout << "marrow.mhd written" << std::endl;
		*/
		
	}

	void ImageHandlerFloat::ITKtoVTKtoITK(){
		
		typedef itk::MinimumMaximumImageCalculator< ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_image);
		calculator->Compute();
		std::cout << "min intensity in ITKtoVTKtoITK - ITK image input: " << calculator->GetMinimum() << std::endl;
		std::cout << "max intensity in ITKtoVTKtoITK - ITK image input: " << calculator->GetMaximum() << std::endl;
	

		// ITK to VTK 
		typedef itk::ImageToVTKImageFilter< ImageType > imageToVTKImageFilter;
		imageToVTKImageFilter::Pointer filter = imageToVTKImageFilter::New();
		filter->SetInput(_image);
		filter->Update();
		_vtkImage->DeepCopy(filter->GetOutput());
		_vtkImage->Update();

		// VTK to ITK
		typedef itk::VTKImageToImageFilter< ImageType > VTKToImageImageFilter;
		VTKToImageImageFilter::Pointer filterBack = VTKToImageImageFilter::New();
		filterBack->SetInput(_vtkImage);
		filterBack->Update();
		ImageType::ConstPointer itk_image = filterBack->GetOutput(); // the filter returns a ConstPointer
		
		calculator->SetImage(itk_image);
		calculator->Compute();
		//std::cout << "min intensity in ITKtoVTKtoITK - ITK image output: " << calculator->GetMinimum() << std::endl;
		//std::cout << "max intensity in ITKtoVTKtoITK - ITK image output: " << calculator->GetMaximum() << std::endl;
	
		typedef itk::ImageRegionIterator< ImageType> IteratorType;
		IteratorType it1 (_image, _image->GetRequestedRegion());
		typedef itk::ImageRegionConstIterator< ImageType> ConstIteratorType;
		ConstIteratorType it2 (itk_image, itk_image->GetRequestedRegion());

		for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2){
			it1.Set(it2.Get());
		}
		
		// write ITK image
		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( _imageFileName.toAscii().data() );
		writer->SetInput( _image );
		//writer->Update();
		std::cout << "metafile written" << std::endl;
	}

	void ImageHandlerFloat::Extrusion(){

		std::cout << std::endl;
		std::cout << "-- Extrusion --" << std::endl;
		
		typedef itk::MinimumMaximumImageCalculator< ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_image);
		calculator->Compute();
		//std::cout << "min intensity in Extrusion: " << calculator->GetMinimum() << std::endl;
		//std::cout << "max intensity in Extrusion: " << calculator->GetMaximum() << std::endl;
	



		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( "beforeExtrusion.mhd" );
		writer->SetInput( _image );
		//writer->Update();


		/****************** ITERATORS ******************/

		// neighboor iterator (connection 26)
		typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
		NeighborhoodIteratorType::RadiusType radius;
		radius.Fill(1);
		
		NeighborhoodIteratorType::OffsetType voxelPosition = {{0,0,0}};
		NeighborhoodIteratorType::OffsetType offset1 = {{-1,-1,0}};
		NeighborhoodIteratorType::OffsetType offset2 = {{1,-1,0}};
		NeighborhoodIteratorType::OffsetType offset3 = {{-1,0,0}};
		NeighborhoodIteratorType::OffsetType offset4 = {{1,0,0}};
		NeighborhoodIteratorType::OffsetType offset5 = {{-1,1,0}};
		NeighborhoodIteratorType::OffsetType offset6 = {{1,1,0}};
		NeighborhoodIteratorType::OffsetType offset7 = {{0,-1,0}};
		NeighborhoodIteratorType::OffsetType offset8 = {{0,1,0}};

		NeighborhoodIteratorType::OffsetType offset9 = {{0,0,-1}};
		NeighborhoodIteratorType::OffsetType offset10 = {{-1,-1,-1}};
		NeighborhoodIteratorType::OffsetType offset11 = {{1,-1,-1}};
		NeighborhoodIteratorType::OffsetType offset12 = {{-1,0,-1}};
		NeighborhoodIteratorType::OffsetType offset13 = {{1,0,-1}};
		NeighborhoodIteratorType::OffsetType offset14 = {{-1,1,-1}};
		NeighborhoodIteratorType::OffsetType offset15 = {{1,1,-1}};
		NeighborhoodIteratorType::OffsetType offset16 = {{0,-1,-1}};
		NeighborhoodIteratorType::OffsetType offset17 = {{0,1,-1}};

		NeighborhoodIteratorType::OffsetType offset18 = {{0,0,1}};
		NeighborhoodIteratorType::OffsetType offset19 = {{-1,-1,1}};
		NeighborhoodIteratorType::OffsetType offset20 = {{1,-1,1}};
		NeighborhoodIteratorType::OffsetType offset21 = {{-1,0,1}};
		NeighborhoodIteratorType::OffsetType offset22 = {{1,0,1}};
		NeighborhoodIteratorType::OffsetType offset23 = {{-1,1,1}};
		NeighborhoodIteratorType::OffsetType offset24 = {{1,1,1}};
		NeighborhoodIteratorType::OffsetType offset25 = {{0,-1,1}};
		NeighborhoodIteratorType::OffsetType offset26 = {{0,1,1}};

		// region iterator
		typedef itk::ImageRegionIterator< ImageType> IteratorType;
		
				
		/****************** MAKING BACKGROUND UNIFORM ******************/
		// expecially for recreated images, the background could be not exactly == -1024. 
		// background voxel value (the voxel in position 2,2,2, is considered having the background value)
		ImageType::IndexType pixelIndex;
		pixelIndex[0] = 15; pixelIndex[1] = 15; pixelIndex[2] = 15;
		ImageType::PixelType backgroundValue = _image->GetPixel( pixelIndex );
		std::cout << "backgroundValue: " << backgroundValue << std::endl;
		backgroundValue = -1024; // it has to be anyway lower than the marrow
		std::cout << "imposed backgroundValue: " << backgroundValue << std::endl;
		/*
		IteratorType itTemp (_image, _image->GetRequestedRegion());
		for (itTemp.GoToBegin(); !itTemp.IsAtEnd(); ++itTemp){
			if (itTemp.Get() < -1000)
				itTemp.Set(backgroundValue); 
			//if (itTemp.Get() > (backgroundValue+100) && itTemp.Get() < (backgroundValue-100))
			//	itTemp.Set(backgroundValue);
		}
		*/
		
		
		/****************** NEEDED IMAGES ******************/
		// copies are needed in order to avoid to take into account what has been changed in one loop during the same current loop

		// original image backup
		ImageType::Pointer originalImage = ImageType::New();
		originalImage->SetRegions( _image->GetRequestedRegion() );
		originalImage->CopyInformation( _image );
		originalImage->Allocate();
		
		// peeled image allocation 
		ImageType::Pointer peeledImage = ImageType::New();
		peeledImage->SetRegions( _image->GetRequestedRegion() );
		peeledImage->CopyInformation( _image );
		peeledImage->Allocate();
	
		// extruded image allocation
		ImageType::Pointer extrudedImage = ImageType::New();
		extrudedImage->SetRegions( _image->GetRequestedRegion() );
		extrudedImage->CopyInformation( _image );
		extrudedImage->Allocate();


		// images iterators
		NeighborhoodIteratorType it1( radius, _image, _image->GetRequestedRegion() );// original image iterator
		NeighborhoodIteratorType it2a( radius, peeledImage, peeledImage->GetRequestedRegion() );// peeled image iterator (the neighboor one is needed)
			
		IteratorType it0 (originalImage, originalImage->GetRequestedRegion());
		IteratorType it1a (_image, _image->GetRequestedRegion());
		IteratorType it2 (peeledImage, peeledImage->GetRequestedRegion());
		IteratorType it3 (extrudedImage, extrudedImage->GetRequestedRegion());
		

		// it keeps the original image in order to compare if the extruded gray level is greater than the original one
		for (it0.GoToBegin(), it1.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1){
			it0.Set( it1.GetPixel(voxelPosition));
		}

		// giving the same background value to the image to peel
		for (it1a.GoToBegin(); !it1a.IsAtEnd(); ++it1a){
			
			if (it1a.Get() < backgroundValue)
				it1a.Set(backgroundValue);
		
		}



		/****************** PEELING ******************/
		// using it1, it1a and it2 -> _image and peeledImage
			
		// copy original image to the peeled one (a simple assignment doesn't allow the two images to be indipendent. changes on one affect also the other)
		for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2){
			it2.Set( it1.GetPixel(voxelPosition));
		}
				
		// peeling off the voxels that are close to the background ones
		ImageType::PixelType value1; ImageType::PixelType value2; ImageType::PixelType value3; ImageType::PixelType value4;
		ImageType::PixelType value5; ImageType::PixelType value6; ImageType::PixelType value7; ImageType::PixelType value8;
		ImageType::PixelType value9; ImageType::PixelType value10; 
		ImageType::PixelType value11; ImageType::PixelType value12; ImageType::PixelType value13; ImageType::PixelType value14;
		ImageType::PixelType value15; ImageType::PixelType value16; ImageType::PixelType value17; ImageType::PixelType value18;
		ImageType::PixelType value19; ImageType::PixelType value20; 
		ImageType::PixelType value21; ImageType::PixelType value22; ImageType::PixelType value23; ImageType::PixelType value24;
		ImageType::PixelType value25; ImageType::PixelType value26; 
		
		// peeling off i times
//////////////////////////////////////
		int layers = 3;
		for (int i=0; i<layers; i++){ 
		
			for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2){

				ImageType::PixelType currentValue;
				currentValue = it1.GetPixel(voxelPosition);
				
				if (currentValue != backgroundValue){
					
					value1 = it1.GetPixel(offset1); value2 = it1.GetPixel(offset2); value3 = it1.GetPixel(offset3); value4 = it1.GetPixel(offset4); 
					value5 = it1.GetPixel(offset5); value6 = it1.GetPixel(offset6); value7 = it1.GetPixel(offset7); value8 = it1.GetPixel(offset8); 
					value9 = it1.GetPixel(offset9); value10 = it1.GetPixel(offset10); 
					value11 = it1.GetPixel(offset11); value12 = it1.GetPixel(offset12); value13 = it1.GetPixel(offset13); value14 = it1.GetPixel(offset14); 
					value15 = it1.GetPixel(offset15); value16 = it1.GetPixel(offset16); value17 = it1.GetPixel(offset17); value18 = it1.GetPixel(offset18); 
					value19 = it1.GetPixel(offset19); value20 = it1.GetPixel(offset20); 
					value21 = it1.GetPixel(offset21); value22 = it1.GetPixel(offset22); value23 = it1.GetPixel(offset23); value24 = it1.GetPixel(offset24); 
					value25 = it1.GetPixel(offset25); value26 = it1.GetPixel(offset26); 
					
					if (value1 == backgroundValue || value2 == backgroundValue || value3 == backgroundValue || value4 == backgroundValue ||
						value5 == backgroundValue || value6 == backgroundValue || value7 == backgroundValue || value8 == backgroundValue ||
						value9 == backgroundValue || value10 == backgroundValue ||
						value11 == backgroundValue || value12 == backgroundValue || value13 == backgroundValue || value14 == backgroundValue ||
						value15 == backgroundValue || value16 == backgroundValue || value17 == backgroundValue || value18 == backgroundValue ||
						value19 == backgroundValue || value20 == backgroundValue ||
						value21 == backgroundValue || value22 == backgroundValue || value23 == backgroundValue || value24 == backgroundValue ||
						value25 == backgroundValue || value26 == backgroundValue){

						it2.Set(backgroundValue);

					}
				}
			}

			// reassigning the peeled image to the original one in order to continue the loop
			for (it1a.GoToBegin(), it2.GoToBegin(); !it1a.IsAtEnd(); ++it1a, ++it2){
				it1a.Set( it2.Get());
			}

			i++; // layers counter
		}

		
		writer->SetFileName( "peeled.mhd" );
		writer->SetInput( peeledImage );
		//writer->Update();


		
		/****************** EXTRUSION ******************/
		// using it0, it2a and it3 
				
		// copy peeled image in the extruded one (a simple assignment doesn't allow the two images to be indipendent. changes on one affect also the other)
		for (it2.GoToBegin(), it3.GoToBegin(); !it2.IsAtEnd(); ++it2, ++it3){
			it3.Set( it2.Get());
		}

		// surrounding with expanded border
////////////////////////
		int extraRing = 2; 
		for (int i=0; i<layers+extraRing; i++){ 

			for (it0.GoToBegin(), it2a.GoToBegin(), it3.GoToBegin(); !it2a.IsAtEnd(); ++it0, ++it2a, ++it3){
			
				ImageType::PixelType currentValue;
				currentValue = it2a.GetPixel(voxelPosition);
								
				if (currentValue == backgroundValue){

					value1 = it2a.GetPixel(offset1); value2 = it2a.GetPixel(offset2); value3 = it2a.GetPixel(offset3); value4 = it2a.GetPixel(offset4); 
					value5 = it2a.GetPixel(offset5); value6 = it2a.GetPixel(offset6); value7 = it2a.GetPixel(offset7); value8 = it2a.GetPixel(offset8); 
					value9 = it2a.GetPixel(offset9); value10 = it2a.GetPixel(offset10); 
					value11 = it2a.GetPixel(offset11); value12 = it2a.GetPixel(offset12); value13 = it2a.GetPixel(offset13); value14 = it2a.GetPixel(offset14); 
					value15 = it2a.GetPixel(offset15); value16 = it2a.GetPixel(offset16); value17 = it2a.GetPixel(offset17); value18 = it2a.GetPixel(offset18); 
					value19 = it2a.GetPixel(offset19); value20 = it2a.GetPixel(offset20); 
					value21 = it2a.GetPixel(offset21); value22 = it2a.GetPixel(offset22); value23 = it2a.GetPixel(offset23); value24 = it2a.GetPixel(offset24); 
					value25 = it2a.GetPixel(offset25); value26 = it2a.GetPixel(offset26); 
					
					if (value1 != backgroundValue || value2 != backgroundValue || value3 != backgroundValue || value4 != backgroundValue ||
						value5 != backgroundValue || value6 != backgroundValue || value7 != backgroundValue || value8 != backgroundValue ||
						value9 != backgroundValue || value10 != backgroundValue ||
						value11 != backgroundValue || value12 != backgroundValue || value13 != backgroundValue || value14 != backgroundValue ||
						value15 != backgroundValue || value16 != backgroundValue || value17 != backgroundValue || value18 != backgroundValue ||
						value19 != backgroundValue || value20 != backgroundValue ||
						value21 != backgroundValue || value22 != backgroundValue || value23 != backgroundValue || value24 != backgroundValue ||
						value25 != backgroundValue || value26 != backgroundValue){

							/*
							// get nonzerovalues and average them
							double value = 0.0; int count = 0.0;
							if (value1 != backgroundValue){value += value1; count++; }// std::cout << value1 << std::endl;}
							if (value2 != backgroundValue){value += value2; count++; }// std::cout << value2 << std::endl;}
							if (value3 != backgroundValue){value += value3; count++; }// std::cout << value3 << std::endl;}
							if (value4 != backgroundValue){value += value4; count++; }// std::cout << value4 << std::endl;}
							if (value5 != backgroundValue){value += value5; count++; }// std::cout << value5 << std::endl;}
							if (value6 != backgroundValue){value += value6; count++; }// std::cout << value6 << std::endl;}
							if (value7 != backgroundValue){value += value7; count++; }// std::cout << value7 << std::endl;}
							if (value8 != backgroundValue){value += value8; count++; }// std::cout << value8 << std::endl;}
							if (value9 != backgroundValue){value += value9; count++; }// std::cout << value9 << std::endl;}
							if (value10 != backgroundValue){value += value10; count++; }// std::cout << value10 << std::endl;}
							if (value11 != backgroundValue){value += value11; count++; }// std::cout << value11 << std::endl;}
							if (value12 != backgroundValue){value += value12; count++; }// std::cout << value12 << std::endl;}
							if (value13 != backgroundValue){value += value13; count++; }// std::cout << value13 << std::endl;}
							if (value14 != backgroundValue){value += value14; count++; }// std::cout << value14 << std::endl;}
							if (value15 != backgroundValue){value += value15; count++; }// std::cout << value15 << std::endl;}
							if (value16 != backgroundValue){value += value16; count++; }// std::cout << value16 << std::endl;}
							if (value17 != backgroundValue){value += value17; count++; }// std::cout << value17 << std::endl;}
							if (value18 != backgroundValue){value += value18; count++; }// std::cout << value18 << std::endl;}
							if (value19 != backgroundValue){value += value19; count++; }// std::cout << value19 << std::endl;}
							if (value20 != backgroundValue){value += value20; count++; }// std::cout << value20 << std::endl;}
							if (value21 != backgroundValue){value += value21; count++; }// std::cout << value21 << std::endl;}
							if (value22 != backgroundValue){value += value22; count++; }// std::cout << value22 << std::endl;}
							if (value23 != backgroundValue){value += value23; count++; }// std::cout << value23 << std::endl;}
							if (value24 != backgroundValue){value += value24; count++; }// std::cout << value24 << std::endl;}
							if (value25 != backgroundValue){value += value25; count++; }// std::cout << value25 << std::endl;}
							if (value26 != backgroundValue){value += value26; count++; }// std::cout << value26 << std::endl;}
							
							value /= count;
							*/

							// take the max value
							float value = backgroundValue;
							if (value1 > value) {value = value1;}
							if (value2 > value) {value = value2;}
							if (value3 > value) {value = value3;}
							if (value4 > value) {value = value4;}
							if (value5 > value) {value = value5;}
							if (value6 > value) {value = value6;}
							if (value7 > value) {value = value7;}
							if (value8 > value) {value = value8;}
							if (value9 > value) {value = value9;}
							if (value10 > value) {value = value10;}
							if (value11 > value) {value = value11;}
							if (value12 > value) {value = value12;}
							if (value13 > value) {value = value13;}
							if (value14 > value) {value = value14;}
							if (value15 > value) {value = value15;}
							if (value16 > value) {value = value16;}
							if (value17 > value) {value = value17;}
							if (value18 > value) {value = value18;}
							if (value19 > value) {value = value19;}
							if (value20 > value) {value = value20;}
							if (value21 > value) {value = value21;}
							if (value22 > value) {value = value22;}
							if (value23 > value) {value = value23;}
							if (value24 > value) {value = value24;}
							if (value25 > value) {value = value25;}
							if (value26 > value) {value = value26;}

							
							// checking if the new value is greater than the original one
							if (value > it0.Get())
								it3.Set(value);
							else
								it3.Set(it0.Get());
					}
				}
			}
			
			// reassing the extruded image to the peeled one in order to continue the loop
			for (it2.GoToBegin(), it3.GoToBegin(); !it2.IsAtEnd(); ++it2, ++it3){
				it2.Set( it3.Get());
			}

			i++; // layers counter
		}
			
		// to output
		_image = extrudedImage;

		writer->SetFileName( "extrudedIn.mhd" );
		writer->SetInput( extrudedImage );
		//writer->Update();


		
	}








	void ImageHandlerFloat::RotateImage(){

		/*
		// getting points
		double center[3]; double pt1[3]; double pt2[3];
		_points->GetPoint(2, pt0);
		_points->GetPoint(0, pt1);
		_points->GetPoint(1, pt2);
		// points to the origin (xAxis = pt2-center; yAxis = pt1-center)
		double pt1origin[3];  pt1origin[0]=pt1[0]-pt0[0]; pt1origin[1]=pt1[1]-pt0[1]; pt1origin[2]=pt1[2]-pt0[2];
		double pt2origin[3];  pt2origin[0]=pt2[0]-pt0[0]; pt2origin[1]=pt2[1]-pt0[1]; pt2origin[2]=pt2[2]-pt0[2];
		// calculation of the third point (vector product)
		double pt3origin[3]; 
		pt3origin[0] = pt1origin[1]*pt2origin[2] - pt1origin[2]*pt2origin[1];
		pt3origin[1] = pt1origin[2]*pt2origin[0] - pt1origin[0]*pt2origin[2];
		pt3origin[2] = pt1origin[0]*pt2origin[1] - pt1origin[1]*pt2origin[0];
		
		// rotation axis norms
		double normPt1 = vtkMath::Norm(pt1origin);
		double normPt2 = vtkMath::Norm(pt2origin);
		double normPt3 = vtkMath::Norm(pt3origin);
		
		// rotation versors
		for(unsigned int i = 0; i < 3; ++i) {
			pt1origin[i]=pt1origin[i]/normPt1;
			pt2origin[i]=pt2origin[i]/normPt2;
			pt3origin[i]=pt3origin[i]/normPt3;
		}	



		// translate image to the origin
		const ImageType::PointType& origin = _image->GetOrigin();
		
		// rotation
		AffineTransformType::Pointer transform = AffineTransformType::New();
		AffineTransformType::MatrixType matrix;
		matrix[0][0] = pt1origin[0]; matrix[0][1] = 0; matrix[0][2] = 0;
		matrix[1][0] = pt1origin[1]; matrix[1][1] = 0; matrix[1][2] = 0;
		matrix[2][0] = pt1origin[2]; matrix[2][1] = 0; matrix[2][2] = 0;
		// center of rotation
		AffineTransformType::InputPointType center;
		center[0] = pt0[0];
		center[1] = pt0[1];
		center[2] = pt0[2];
		// translation
		AffineTransformType::OutputVectorType translate;
		translate[0] = origin[0];
		translate[1] = origin[1];
		translate[2] = origin[2];
		
		transform->SetMatrix(matrix);
		transform->SetCenter(center);
		transform->SetTranslation(translate);


		// resample label image
		typedef itk::ResampleImageFilter< ImageType, ImageType> ResampleFilterType;
		ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
		resampleFilter->SetInput(m_MaskImage);
		resampleFilter->SetTransform(transform);
		resampleFilter->UseReferenceImageOn();
		resampleFilter->SetReferenceImage(m_MaskImage);
		resampleFilter->SetDefaultPixelValue(AbstractRegistration::OUTSIDE_LABEL);
		
		typedef itk::NearestNeighborInterpolateImageFunction< AbstractRegistration::LabelImageType > LabelInterpolatorType;
		LabelInterpolatorType::Pointer labelInterpolator = LabelInterpolatorType::New();
		labelResampler->SetInterpolator(labelInterpolator);
		labelResampler->Update();
		
		*/

	}
	void ImageHandlerFloat::ExtractVolumeImage(){

		////////////////////////////////////////////////////////////
		// allocate a third image (neck intensities) for test saving
		/*ImageType::Pointer test = ImageType::New();
		test->SetRegions( _image->GetRequestedRegion() );
		test->CopyInformation( _image );
		test->Allocate();
		test->FillBuffer(-1024);
		*////////////////////////////////////////////////////////////


		// iterator for the mask
		typedef itk::ImageRegionIteratorWithIndex< ImageType> IteratorType;
		IteratorType it (_mask, _mask->GetRequestedRegion());
		
		// vnl vectors want allocation -> vectorTemp first
		std::vector<double> vectorTemp;
		
		for (it.GoToBegin(); !it.IsAtEnd(); ++it){

			if (it.Get() == 1){
				
				// get the white voxel index
				ImageType::IndexType maskIndex;
				maskIndex = it.GetIndex();
				//std::cout << "maskIndex: " << maskIndex[0] << ' ' << maskIndex[1] << ' ' << maskIndex[2] << std::endl;   
				
				// look in the svf
				FieldType::PixelType svfValue = _svf->GetPixel( maskIndex );
				//std::cout << "svfValue : " << svfValue [0] << ' ' << svfValue [1] << ' ' << svfValue [2] << std::endl;

				// calculate index for the calibrated image
				ImageType::IndexType imageIndex;
				imageIndex[0] = floor (maskIndex[0] + svfValue[0] + 0.5);
				imageIndex[1] = floor (maskIndex[1] + svfValue[1] + 0.5);
				imageIndex[2] = floor (maskIndex[2] + svfValue[2] + 0.5);
				//std::cout << "imageIndex: " << imageIndex[0] << ' ' << imageIndex[1] << ' ' << imageIndex[2] << std::endl;


				// extract the value in the calibrated image
				ImageType::PixelType imageValue = _image->GetPixel( imageIndex );
				//std::cout << "imageValue : " << imageValue << std::endl;

				vectorTemp.push_back(imageValue);
								
				
				////////////////////////////////////////////////////////////
				//test->SetPixel( imageIndex, imageValue ); 
				////////////////////////////////////////////////////////////


			}
		}


		// vnl vector
		_vector.set_size(vectorTemp.size());
		for (int i=0; i<vectorTemp.size(); i++){
			_vector(i) = vectorTemp.at(i);
		}		


		////////////////////////////////////////////////////////////
		/*typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( "test.mhd" );
		writer->SetInput( test );
		writer->Update();
		std::cout << "metafile written" << std::endl;
		*////////////////////////////////////////////////////////////

	}
}