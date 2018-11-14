/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <ImageHandler.h>

#include <algorithm>

#include <StatisticsEntropy.h>
#include <VectorImageHandler.h>
#include <VnlWriterMatrix.h>
#include <VnlWriterVector.h>

#include <itkAddConstantToImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkHistogram.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageSeriesReader.h>
#include <itkImageSeriesWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToVTKImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkListSampleToHistogramGenerator.h>
#include <itkMedianImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkScalarImageToListAdaptor.h>
#include <itkVTKImageToImageFilter.h>
#include <itkWarpImageFilter.h>

#include <vtkDoubleArray.h>
#include <vtkMarchingCubes.h>

#include <vnl/vnl_vector.h>

#include <QFile>
#include <QTextStream>

using namespace statistics;
using namespace vnl;

namespace image {

	// constructor
	ImageHandler::ImageHandler(){

		_spacing[0] = 1.25; _spacing[1] = 1.25; _spacing[2] = 1.25;
		_threshold = -1000;
		_edgePaddingValue = -1024;

		_mesh = vtkPolyData::New();
		_vtkImage = vtkImageData::New();
			
	}

	// destructor
	ImageHandler::~ImageHandler(){
		
		_mesh->Delete();
		_vtkImage->Delete();
		_image = NULL;


	}

	
	// member functions
	
	/******************************** IMAGE READERS ********************************/
	void ImageHandler::DicomReaderUpdate(){

		// uncomprehensible commands from itk
		typedef itk::ImageSeriesReader< ImageHandler::ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		
		typedef itk::GDCMImageIO ImageIOType;
		ImageIOType::Pointer dicomIO = ImageIOType::New();
		reader->SetImageIO( dicomIO );

		typedef itk::GDCMSeriesFileNames NamesGeneratorType;
		NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
		nameGenerator->SetUseSeriesDetails( true );
		nameGenerator->SetDirectory( _imageFolderFileName.ascii() );

		typedef std::vector< std::string > SeriesIdContainer;
		const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
		SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
		SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
		while( seriesItr != seriesEnd )
		{
		//std::cout << seriesItr->c_str() << std::endl;
		seriesItr++;
		}
		std::string seriesIdentifier;
		seriesIdentifier = seriesUID.begin()->c_str();

		typedef std::vector< std::string > FileNamesContainer;
		FileNamesContainer fileNames;
		fileNames = nameGenerator->GetFileNames( seriesIdentifier );
		typedef itk::ImageSeriesReader< ImageHandler::ImageType > ImageSeriesReaderType;
		ImageSeriesReaderType::Pointer imageSeriesReader = ImageSeriesReaderType::New();
		imageSeriesReader->SetFileNames( fileNames );
		imageSeriesReader->Update();

		// image
		_image = imageSeriesReader->GetOutput();	
		
		// image data
		const ImageHandler::ImageType::SpacingType& spacing = _image->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = _image->GetOrigin();
		const ImageHandler::ImageType::SizeType& size = _image->GetBufferedRegion().GetSize();
		std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
		std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
		std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
		typedef itk::MinimumMaximumImageCalculator< ImageHandler::ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_image);
		calculator->Compute();
		std::cout << "min intensity: " << calculator->GetMinimum() << "; max intensity: " << calculator->GetMaximum() << std::endl;
		
	}

	void ImageHandler::DicomMaskReaderUpdate(){

		// uncomprehensible commands from itk
		typedef itk::ImageSeriesReader< ImageHandler::ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		
		typedef itk::GDCMImageIO ImageIOType;
		ImageIOType::Pointer dicomIO = ImageIOType::New();
		reader->SetImageIO( dicomIO );

		typedef itk::GDCMSeriesFileNames NamesGeneratorType;
		NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
		nameGenerator->SetUseSeriesDetails( true );
		nameGenerator->SetDirectory( _maskFolderFileName.ascii() );

		typedef std::vector< std::string > SeriesIdContainer;
		const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
		SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
		SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

		while( seriesItr != seriesEnd )
		{
		//std::cout << seriesItr->c_str() << std::endl;
		seriesItr++;
		}
		std::string seriesIdentifier;
		seriesIdentifier = seriesUID.begin()->c_str();

		typedef std::vector< std::string > FileNamesContainer;
		FileNamesContainer fileNames;
		fileNames = nameGenerator->GetFileNames( seriesIdentifier );
		typedef itk::ImageSeriesReader< ImageHandler::ImageType > ImageSeriesReaderType;
		ImageSeriesReaderType::Pointer imageSeriesReader = ImageSeriesReaderType::New();
		imageSeriesReader->SetFileNames( fileNames );
		imageSeriesReader->Update();

		// mask
		_mask = imageSeriesReader->GetOutput();	

		// mask data
		const ImageHandler::ImageType::SpacingType& spacing = _mask->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = _mask->GetOrigin();
		const ImageHandler::ImageType::SizeType& size = _mask->GetBufferedRegion().GetSize();
		std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
		std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
		std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
		typedef itk::MinimumMaximumImageCalculator< ImageHandler::ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_mask);
		calculator->Compute();
		std::cout << "min intensity: " << calculator->GetMinimum() << "; max intensity: " << calculator->GetMaximum() << std::endl; 
		

	}
	void ImageHandler::MetafileReaderUpdate(){
				
		typedef itk::ImageFileReader< ImageHandler::ImageType > ReaderType;
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
		const ImageHandler::ImageType::SpacingType& spacing = _image->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = _image->GetOrigin();
		const ImageHandler::ImageType::SizeType& size = _image->GetBufferedRegion().GetSize();
		std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
		std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
		std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
		typedef itk::MinimumMaximumImageCalculator< ImageHandler::ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_image);
		calculator->Compute();
		std::cout << "min intensity: " << calculator->GetMinimum() << "; max intensity: " << calculator->GetMaximum() << std::endl; 
	
	}

	void ImageHandler::MaskReaderUpdate() {

		typedef itk::ImageFileReader< ImageHandler::ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(_maskFileName.toAscii().data());
		//std::cout << _imageFileName.toAscii().data() << std::endl;
		try{
			reader->Update();
		}
		catch( itk::ExceptionObject & excp ){
			std::cerr << "!!! Problem reading the input file" << std::endl;
		}
		
		// mask
		_mask = reader->GetOutput();
		
		// mask data
		const ImageHandler::ImageType::SpacingType& spacing = _image->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = _image->GetOrigin();
		const ImageHandler::ImageType::SizeType& size = _image->GetBufferedRegion().GetSize();
		std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
		std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
		std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
		typedef itk::MinimumMaximumImageCalculator< ImageHandler::ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_mask);
		calculator->Compute();
		std::cout << "min intensity: " << calculator->GetMinimum() << "; max intensity: " << calculator->GetMaximum() << std::endl; 
	
	}

	void ImageHandler::FieldReaderUpdate(){

		typedef itk::ImageFileReader< ImageHandler::FieldType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(_fieldFileName.toAscii().data());
		//std::cout << _fieldFileName.toAscii().data() << std::endl;
		try{
			reader->Update();
		}
		catch( itk::ExceptionObject & excp ){
			std::cerr << "!!! Problem reading the input file" << std::endl;
		}
		
		// image
		_field = reader->GetOutput();
		
		// image data
		const ImageHandler::ImageType::SpacingType& spacing = _field->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = _field->GetOrigin();
		const ImageHandler::ImageType::SizeType& size = _field->GetBufferedRegion().GetSize();
		std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
		std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
		std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
		
	}

	/******************************** IMAGE WRITERS ********************************/
	void ImageHandler::MetafileWriterUpdate(){

		typedef itk::ImageFileWriter< ImageHandler::ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( _imageFileName.toAscii().data() );
		writer->SetInput( _image );
		writer->Update();
		std::cout << "metafile written" << std::endl;
	
	}

	void ImageHandler::DicomWriterUpdate(){

		// typedefs
		typedef itk::GDCMImageIO ImageIOType;
		ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
		typedef signed short OutputPixelType;
		const unsigned int OutputDimension = 2;
		typedef itk::Image< OutputPixelType, OutputDimension > Image2DType;
		typedef itk::ImageSeriesWriter<	ImageHandler::ImageType, Image2DType > SeriesWriterType;
		SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
		seriesWriter->SetImageIO( gdcmImageIO );
		ImageHandler::ImageType::SizeType inputSize2 = _image->GetBufferedRegion().GetSize();
		std::vector<std::string> fileNames;

		QString path; 
		for(int a=0; a<inputSize2[2]; a++){ 
			path = _imageFolderFileName;
			if ((a+1)<10)
				path.append("/00");
			else if ((a+1)>=10 && (a+1)<100)
				path.append("/0");
			else
				path.append("/");
			QString temp= QString("%1").arg(a+1);
			path.append(temp);
			path.append(".dcm");
			fileNames.push_back(path.toAscii().data());
			
			//std::cout << path.toAscii().data() << std::endl;

		}
		
		/*
		ImageHandler::ImageType::IndexType pixelIndex;
		pixelIndex[0] = 1; // x position
		pixelIndex[1] = 1; // y position
		pixelIndex[2] = 1; // z position

		ImageHandler::ImageType::PixelType pixelValue = _image->GetPixel( pixelIndex );
		std::cout << pixelValue << std::endl;
		*/

		seriesWriter->SetInput( _image );
		seriesWriter->SetFileNames( fileNames );
		try{
			seriesWriter->Update();
		}
		catch( itk::ExceptionObject & excp ){
		std::cerr << "Exception thrown while writing the series " << std::endl;
		std::cerr << excp << std::endl;
		}
	}
	
	/******************************** IMAGE MASKING ********************************/
	void ImageHandler::MaskToImage(){
		
		ImageHandler::ImageType::IndexType pixelIndex;
		pixelIndex[0] = 1; // x position
		pixelIndex[1] = 1; // y position
		pixelIndex[2] = 1; // z position
		ImageHandler::ImageType::PixelType pixelValue = _image->GetPixel( pixelIndex );
		
		typedef itk::ImageRegionIterator< ImageHandler::ImageType> IteratorType;	
		IteratorType it1 (_image, _image->GetRequestedRegion());
		IteratorType it2 (_mask, _mask->GetRequestedRegion());
		
		for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2){

			if (it2.Get() == 0){ // mask background
				//it1.Set(pixelValue);
				it1.Set(-1024);
			}
		}
	}
	void ImageHandler::MaskGreyLevelsUpdate(){
		
		// neighboor filter
		typedef itk::ConstNeighborhoodIterator< ImageHandler::ImageType > NeighborhoodIteratorType;
		NeighborhoodIteratorType::RadiusType radius;
		radius.Fill(1);
		NeighborhoodIteratorType it( radius, _image, _image->GetRequestedRegion() );

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


		// give a value that is 100, 200, 300 to the transition voxels based on its neighboothood
		typedef itk::ImageRegionIterator< ImageHandler::ImageType> IteratorType;	
		IteratorType it2 (_image, _image->GetRequestedRegion());

		for (it.GoToBegin(), it2.GoToBegin(); !it.IsAtEnd(); ++it, ++it2){

			ImageType::PixelType value;
			value = it.GetPixel(voxelPosition);
			
			if (value != 300 && value != 200 && value != 100 && value != 0){
				//std::cout << "value: " << value << std::endl;
				
				int c = 0; int t = 0; int m = 0;
				
				// same slide
				value = it.GetPixel(offset1);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset2);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset3);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset4);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset5);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset6);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset7);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset8);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				
				// slide below
				value = it.GetPixel(offset9);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset10);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset11);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset12);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset13);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset14);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset15);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset16);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset17);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				
				// slide above
				value = it.GetPixel(offset18);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset19);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset20);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset21);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset22);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset23);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset24);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset25);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone
				value = it.GetPixel(offset26);
				if (value == 300) c++; else if (value == 200) t++; else if (value == 100) m++; else if (value == 0) c++; // if part of the background, it will become cortical bone

				// assign the grey level
				if ((c>t && c>m) || (c>=t && c>m) || (c>t && c>=m) || (c>=t && c>=m))
					value = 300;
				if ((t>c && t>m) || (t>=c && t>m) || (t>c && t>=m) || (t>=c && t>=m))
					value = 200;
				if ((m>t && m>c) || (m>=t && m>c) || (m>t && m>=c) || (m>=t && m>=c))
					value = 100;

				it2.Set(value);
			}
		}

		
		// correction of the outliers between 0 and 300 that are between the cortical bone and the background
		ImageType::SizeType size = _image->GetLargestPossibleRegion().GetSize(); 
		
		for (int i=1; i<size[0]-1; i++){
			for (int j=1; j<size[1]-1; j++){
				for (int k=1; k<size[2]-1; k++){

					ImageType::IndexType pixelIndex;
					pixelIndex[0] = i; pixelIndex[1] = j; pixelIndex[2] = k; 

					ImageType::PixelType value;
					value = _image->GetPixel(pixelIndex);

					if (value == 0){
						
						pixelIndex[0] = i-1; pixelIndex[1] = j-1; pixelIndex[2] = k; 
						value = _image->GetPixel(pixelIndex);
						if ((value != 300) && (value != 0)) _image->SetPixel(pixelIndex, 300);

						pixelIndex[0] = i+1; pixelIndex[1] = j-1; pixelIndex[2] = k; 
						value = _image->GetPixel(pixelIndex);
						if ((value != 300) && (value != 0)) _image->SetPixel(pixelIndex, 300);

						pixelIndex[0] = i-1; pixelIndex[1] = j; pixelIndex[2] = k; 
						value = _image->GetPixel(pixelIndex);
						if ((value != 300) && (value != 0)) _image->SetPixel(pixelIndex, 300);

						pixelIndex[0] = i+1; pixelIndex[1] = j; pixelIndex[2] = k; 
						value = _image->GetPixel(pixelIndex);
						if ((value != 300) && (value != 0)) _image->SetPixel(pixelIndex, 300);

						pixelIndex[0] = i-1; pixelIndex[1] = j+1; pixelIndex[2] = k; 
						value = _image->GetPixel(pixelIndex);
						if ((value != 300) && (value != 0)) _image->SetPixel(pixelIndex, 300);

						pixelIndex[0] = i+1; pixelIndex[1] = j+1; pixelIndex[2] = k; 
						value = _image->GetPixel(pixelIndex);
						if ((value != 300) && (value != 0)) _image->SetPixel(pixelIndex, 300);

						pixelIndex[0] = i; pixelIndex[1] = j-1; pixelIndex[2] = k; 
						value = _image->GetPixel(pixelIndex);
						if ((value != 300) && (value != 0)) _image->SetPixel(pixelIndex, 300);

						pixelIndex[0] = i; pixelIndex[1] = j+1; pixelIndex[2] = k; 
						value = _image->GetPixel(pixelIndex);
						if ((value != 300) && (value != 0)) _image->SetPixel(pixelIndex, 300);

					}
				}			
			}
		}
	}
	
	/******************************** ITK VTK CONVERSION ***************************/
	void ImageHandler::ITKtoVTK(){

		typedef itk::ImageToVTKImageFilter< ImageType > imageToVTKImageFilter;
		imageToVTKImageFilter::Pointer filter = imageToVTKImageFilter::New();
		filter->SetInput(_image);
		filter->Update();

		_vtkImage->DeepCopy(filter->GetOutput());
		_vtkImage->Update();

		filter->Delete();
		
	}

	void ImageHandler::ITKtoVTKtoITK(){
		
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
		
		typedef itk::ImageRegionIterator< ImageType> IteratorType;
		IteratorType it1 (_image, _image->GetRequestedRegion());
		typedef itk::ImageRegionConstIterator< ImageType> ConstIteratorType;
		ConstIteratorType it2 (itk_image, itk_image->GetRequestedRegion());

		for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2){
			it1.Set(it2.Get());
		}
		
		// write ITK image
		typedef itk::ImageFileWriter< ImageHandler::ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( _imageFileName.toAscii().data() );
		writer->SetInput( _image );
		//writer->Update();
		std::cout << "metafile written" << std::endl;
	}

	
	/*********************************** PROCESSING ********************************/
	void ImageHandler::MarchingCubesUpdate(){
		
		typedef itk::ImageToVTKImageFilter< ImageType > imageToVTKImageFilter;
		imageToVTKImageFilter::Pointer filter = imageToVTKImageFilter::New();
		filter->SetInput(_image);
		filter->Update();

		_vtkImage->DeepCopy(filter->GetOutput());
		_vtkImage->Update();
						
		vtkMarchingCubes* marchingCubes = vtkMarchingCubes::New();
		marchingCubes->SetInput(_vtkImage);
		marchingCubes->SetValue(0,_threshold); 
		marchingCubes->Update();

		_mesh->DeepCopy(marchingCubes->GetOutput());

		marchingCubes->Delete();

	}


	void ImageHandler::ResampleUpdate(){
	
		// resampler
		typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
		ResampleFilterType::Pointer resampler = ResampleFilterType::New();
		
		// linear interpolator
		typedef itk::LinearInterpolateImageFunction < ImageType, double > InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		resampler->SetInterpolator( interpolator );

		// spacing
		resampler->SetOutputSpacing( _spacing );
		
		// origin
		resampler->SetOutputOrigin( _image->GetOrigin() );
		
		// size
		ImageType::SizeType inputSize = _image->GetLargestPossibleRegion().GetSize(); // input size
		typedef ImageType::SizeType::SizeValueType SizeValueType; // output size
		ImageType::SizeType size;
		for(unsigned int i = 0; i < 3; ++i)
			size[i] = ( inputSize[i] * _image->GetSpacing()[i] ) / _spacing[i];
		resampler->SetSize( size );
	
		// direction
		resampler->SetOutputDirection( _image->GetDirection() );

		// update
		resampler->SetInput(_image);
		resampler->Update();

		ImageType::PointType orig;
		orig = resampler->GetOutput()->GetOrigin();
		std::cout << orig[0] << ' ' << orig[1] << ' ' << orig[2] << std::endl;

		_image = resampler->GetOutput();
				
	}
	void ImageHandler::ApplyDVFToImageUpdate(){

		typedef itk::WarpImageFilter< ImageHandler::ImageType,ImageHandler::ImageType, ImageHandler::FieldType> WarpImageFilterType;
		WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();
		typedef WarpImageFilterType::CoordRepType CoordRepType;
				
		//typedef EdgeLinearInterpolateImageFunction< AbstractRegistration::ImageType,AbstractRegistration::ImageType,CoordRepType > EdgeInterpolatorType;
		//EdgeInterpolatorType::Pointer interpolator = EdgeInterpolatorType::New();
		//interpolator->SetInputMaskImage(_referenceMask);
		//typedef itk::NearestNeighborInterpolateImageFunction< ImageHandler::ImageType, CoordRepType> InterpolatorType;
		//InterpolatorType::Pointer interpolator = InterpolatorType::New();
		typedef itk::LinearInterpolateImageFunction< ImageHandler::ImageType,CoordRepType > InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		
		warpImageFilter->SetInput(_image);
		warpImageFilter->SetDeformationField(_constField);
		warpImageFilter->SetInterpolator(interpolator);
		warpImageFilter->SetOutputSpacing(_image->GetSpacing());
		warpImageFilter->SetOutputOrigin(_image->GetOrigin());
		warpImageFilter->SetEdgePaddingValue(_edgePaddingValue);
		warpImageFilter->Update();
		_image = warpImageFilter->GetOutput();		
	}	
	void ImageHandler::AverageImageGreyLevels(){

		// read first image to get characteristics
		typedef itk::ImageFileReader< ImageHandler::ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(_imageFileNames[0]);
		try{
			reader->Update();
		}
		catch( itk::ExceptionObject & excp ){
			std::cerr << "!!! Problem reading the input file" << std::endl;
		}
		// allocate average image
		ImageHandler::ImageType::Pointer averageImage = ImageType::New();
		averageImage->SetRegions( reader->GetOutput()->GetRequestedRegion() );
		averageImage->CopyInformation( reader->GetOutput() );
		averageImage->Allocate();
		ImageHandler::ImageType::PixelType p = itk::NumericTraits< ImageHandler::ImageType::PixelType >::Zero;
		averageImage->FillBuffer(p);

		const ImageHandler::ImageType::SpacingType& spacing = averageImage->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = averageImage->GetOrigin();
		const ImageHandler::ImageType::SizeType& size = averageImage->GetBufferedRegion().GetSize();
		std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
		std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
		std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;

		typedef itk::AddImageFilter< ImageHandler::ImageType, ImageHandler::ImageType, ImageHandler::ImageType > AddFilterType;
		AddFilterType::Pointer addFilter = AddFilterType::New();

		for (int i=0; i<_imageFileNames.size(); i++){

			// read image
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(_imageFileNames[i]);
			std::cout << "image: " << _imageFileNames[i].toAscii().data() << std::endl;
			try{
				reader->Update();
			}
			catch( itk::ExceptionObject & excp ){
				std::cerr << "!!! Problem reading the input file" << std::endl;
			}

			// read correspondent VF
			VectorImageHandler* vectorHandler = new VectorImageHandler;
			QString temp = _imageFileNames[i];
			temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
			temp.remove(0, temp.lastIndexOf("/")+1);
			QString meshTemp;
			for (int y=0; y<_fieldFileNames.size(); y++){
				if (_fieldFileNames[y].contains(temp)){
					std::cout << "velocity field: " << _fieldFileNames[y].ascii() << std::endl;
					vectorHandler->SetFieldFileName(_fieldFileNames[i].ascii());
					vectorHandler->MetafileReaderUpdate();
				}
			}

			// VF to DVF
			std::cout << "convert VF to DVF" <<  std::endl;
			vectorHandler->VFtoDVFinverted();

			//apply DVF to image
			std::cout << "apply DVF to image" << std::endl;
			typedef itk::WarpImageFilter< ImageHandler::ImageType,ImageHandler::ImageType, ImageHandler::FieldType> WarpImageFilterType;
			WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();
			typedef WarpImageFilterType::CoordRepType CoordRepType;
			typedef itk::LinearInterpolateImageFunction< ImageHandler::ImageType,CoordRepType > InterpolatorType;
			InterpolatorType::Pointer interpolator = InterpolatorType::New();
			warpImageFilter->SetInput(_image);
			warpImageFilter->SetDeformationField(_field);
			warpImageFilter->SetInterpolator(interpolator);
			warpImageFilter->SetOutputSpacing(_image->GetSpacing());
			warpImageFilter->SetOutputOrigin(_image->GetOrigin());
			warpImageFilter->SetEdgePaddingValue(_edgePaddingValue);
			warpImageFilter->Update();

			// save image
			_imageFileNames[i].replace(QString("_velocity_field.mhd"), QString("_warped.mhd"));
			std::cout << "warped image: " << _imageFileNames[i].ascii() << std::endl;
			typedef itk::ImageFileWriter< ImageHandler::ImageType > WriterType;
			WriterType::Pointer writer = WriterType::New();
			writer->SetFileName( _imageFileNames[i].ascii());
			writer->SetInput( warpImageFilter->GetOutput() );
			writer->Update();
			

			// add to average
			typedef itk::MultiplyByConstantImageFilter< ImageHandler::ImageType, double, ImageHandler::ImageType > MultiplyFilter;
			MultiplyFilter::Pointer multiplyFilter = MultiplyFilter::New();
			multiplyFilter->SetInput(warpImageFilter->GetOutput());
			multiplyFilter->SetConstant(1.0 / _imageFileNames.size());
			
			
			addFilter->SetInput1( averageImage );
			addFilter->SetInput2( multiplyFilter->GetOutput() );
		}

		// save average
		QString temp = _imageFileNames[0];
		temp.remove(temp.lastIndexOf("/"),temp.size()-1);
		temp.append("average.mhd");
		std::cout << "average: " << temp.toAscii().data() << std::endl;
		typedef itk::ImageFileWriter< ImageHandler::ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( temp);
		writer->SetInput( addFilter->GetOutput() );
		writer->Update();
}
	/*void ImageHandler::PseudoCalibration(){	

		typedef itk::MinimumMaximumImageCalculator< ImageHandler::ImageType > MinMaxCalculator;
		MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
		calculator->SetImage(_image);
		calculator->Compute();
		
		//double zeroDensity = 0.0; double maxDensity = 1.08; // y-axis // this is rhoAsh
		//double zeroHU = 0.0; double maxHU = calculator->GetMaximum(); // x-axis
		double zeroDensity = 0.15; double maxDensity = 1.08; // y-axis // this is rhoAsh
		double zeroHU = -90.0; double maxHU = calculator->GetMaximum(); // x-axis
		double slope = (maxDensity - zeroDensity)/(maxHU - zeroHU); // slope
		double intercept = zeroDensity - ((zeroDensity-maxDensity)/(zeroHU-maxHU))*zeroHU;
		//double slope = (maxDensity - zeroDensity)/(maxHU - zeroHU); // slope
		//double intercept = 0.0; // intercept
		//double slope = 0.00068105; // Schileo
		//double intercept = 0.0740;
			
		std::cout << "slope: " << slope << "  intercept: " << intercept << std::endl;
		
		typedef itk::ImageRegionIterator< ImageType> IteratorType;
		IteratorType it (_image, _image->GetRequestedRegion());

		ImageType::IndexType pixelIndex;
		pixelIndex[0] = 2; pixelIndex[1] = 2; pixelIndex[2] = 2;
		ImageType::PixelType backgroundValue = _image->GetPixel( pixelIndex );


		for (it.GoToBegin(); !it.IsAtEnd(); ++it){

			double greyLevel = it.Get();
			
			if (greyLevel != backgroundValue){
				greyLevel *= slope;
				greyLevel += intercept;
				greyLevel *= 1000; // multiplication by 1000 in order to keep short images (could be considered as g/mm3)
				it.Set(greyLevel);
			}
		}

	}*/	
	void ImageHandler::Extrusion(){
		
		typedef itk::ImageFileWriter< ImageHandler::ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( "beforeExtrusion.mhd" );
		writer->SetInput( _image );
		//writer->Update();


		/****************** ITERATORS ******************/

		// neighboor iterator (connection 26)
		typedef itk::ConstNeighborhoodIterator< ImageHandler::ImageType > NeighborhoodIteratorType;
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
		backgroundValue = -100;
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
		for (int i=0; i<2; i++){ 
		
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
		writer->Update();


		
		/****************** EXTRUSION ******************/
		// using it2, it2a and it3 
				
		// copy peeled image in the extruded one (a simple assignment doesn't allow the two images to be indipendent. changes on one affect also the other)
		for (it2.GoToBegin(), it3.GoToBegin(); !it2.IsAtEnd(); ++it2, ++it3){
			it3.Set( it2.Get());
		}

		// surrounding with expanded border
		for (int i=0; i<8; i++){ // 3 layers

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
		writer->Update();


		
	}

	void ImageHandler::Crop(){
		
		const ImageHandler::ImageType::SpacingType& spacing = _image->GetSpacing();
		const ImageHandler::ImageType::PointType& origin = _image->GetOrigin();
		ImageHandler::ImageType::SizeType size = _image->GetBufferedRegion().GetSize();

		// input region of interest
		ImageType::RegionType				inputRegion;
		ImageType::RegionType::IndexType	inputStart;
		inputStart[0] = 0; inputStart[1] = 0; inputStart[2] = 0;
		size[2] = 80; // 
		inputRegion.SetSize( size );
		inputRegion.SetIndex( inputStart );
		
		// output region of interest
		ImageType::RegionType outputRegion;
		ImageType::RegionType::IndexType outputStart;
		outputStart[0] = 0;	outputStart[1] = 0;	outputStart[2] = 0;
		size[2] = 100; // 
		outputRegion.SetSize( size );
		outputRegion.SetIndex( outputStart );
		ImageType::Pointer outputImage = ImageType::New();
		outputImage->SetRegions( outputRegion );
		outputImage->SetSpacing( spacing );
		outputImage->SetOrigin( origin );
		outputImage->Allocate();
		ImageType::IndexType pixelIndex;
		pixelIndex[0] = 2; pixelIndex[1] = 2; pixelIndex[2] = 2;
		ImageType::PixelType pixelValue = _image->GetPixel( pixelIndex );
		outputImage->FillBuffer(pixelValue);

		// cropping image
		typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
		typedef itk::ImageRegionIterator< ImageType> IteratorType;
		ConstIteratorType inputIt( _image, inputRegion );
		IteratorType outputIt( outputImage, outputRegion );
		for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt){
			outputIt.Set( inputIt.Get() );
		}

		typedef itk::ImageFileWriter< ImageType > ImageWriterType;
		ImageWriterType::Pointer imageWriter = ImageWriterType::New();
		_imageFileName.replace(QString(".mhd"), QString("_cropped.mhd"));
		imageWriter->SetFileName(_imageFileName);	
		imageWriter->SetInput(outputImage);
		imageWriter->Update();
		
		/*
		// to use in the main if needed
		QString fileName = QFileDialog::getOpenFileName(this, "Load images", " ", "Images (*.mhd)");
		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImageFileName(fileName);
		imageHandler->MetafileReaderUpdate();
		imageHandler->Crop();
		std::cout << "done" << std::endl;
		*/
	}
	void ImageHandler::DicomHeaderReaderFromSlice(){

		// e.g. dicom tags in http://www.sno.phy.queensu.ca/~phil/exiftool/TagNames/DICOM.html
		
		typedef itk::ImageFileReader< ImageType > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		typedef itk::GDCMImageIO ImageIOType;
		ImageIOType::Pointer dicomIO = ImageIOType::New();
		reader->SetFileName(_imageFileName.toAscii().data());
		reader->SetImageIO(dicomIO);
		reader->Update();

		typedef itk::MetaDataDictionary DictionaryType;
		const DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
		typedef itk::MetaDataObject< std::string > MetaDataStringType;

		DictionaryType::ConstIterator end = dictionary.End();

		// rescale slope 
		std::string entryId = "0028|1053";
		DictionaryType::ConstIterator itr = dictionary.Find( entryId );
		std::string tagvalue;
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "Rescale slope: " << tagvalue << std::endl;

		// rescale intercept 
		entryId = "0028|1052";
		itr = dictionary.Find( entryId );
				
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "Rescale intercept: " << tagvalue << std::endl;

		// rescale type 
		entryId = "0028|1054";
		itr = dictionary.Find( entryId );
				
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "Rescale type: " << tagvalue << std::endl;


		// institution name
		entryId = "0008|0080";
				
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "Institution name: " << tagvalue << std::endl;
		
		// manufacturer
		entryId = "0008|0070";
		itr = dictionary.Find( entryId );
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "Manufacturer: " << tagvalue << std::endl;
		
		// manufacturer model name
		entryId = "0008|1090";
		itr = dictionary.Find( entryId );
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "Manufacturer model name: " << tagvalue << std::endl;
		
		// device serial number
		entryId = "0018|1000";
		itr = dictionary.Find( entryId );
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "Device serial number: " << tagvalue << std::endl;	
		
		// KVP
		entryId = "0018|0060";
		itr = dictionary.Find( entryId );
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "KVP: " << tagvalue << std::endl;	
			
		// x-ray current tube
		entryId = "0018|1151";
		itr = dictionary.Find( entryId );
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "x-ray current tube: " << tagvalue << std::endl;

		// pixel spacing
		entryId = "0028|0030";
		itr = dictionary.Find( entryId );
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "pixel spacing: " << tagvalue << std::endl;
		
		// slice thickness
		entryId = "0018|0050";
		itr = dictionary.Find( entryId );
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "slice thickness: " << tagvalue << std::endl;

		// generator power
		entryId = "0018|1170";
		itr = dictionary.Find( entryId );
		
		if( itr != end ){
			MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(itr->second.GetPointer() );
			if( entryvalue ){
				tagvalue = entryvalue->GetMetaDataObjectValue();
			}
		}
		std::cout << "generator power: " << tagvalue << std::endl;
		
		
	}
	void ImageHandler::ImageHistogram(){

		// entropy vector
		vnl_vector<double> entropy;
		entropy.set_size(_imageFileNames.size());
		

		for (int i=0; i<_imageFileNames.size(); i++){
			
			// reading the image
			typedef itk::ImageFileReader< ImageHandler::ImageType > ReaderType;
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(_imageFileNames[i].ascii());
			std::cout <<_imageFileNames[i].ascii() << std::endl;
			try{
				reader->Update();
			}
			catch( itk::ExceptionObject & excp ){
				std::cerr << "!!! Problem reading the input file" << std::endl;
			}
			
			typedef itk::MinimumMaximumImageCalculator< ImageHandler::ImageType > MinMaxCalculator;
			MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
			calculator->SetImage(reader->GetOutput());
			calculator->Compute();
			std::cout << "min intensity: " << calculator->GetMinimum() << "; max intensity: " << calculator->GetMaximum() << std::endl; 
			
			// creating the histogram
			typedef itk::Statistics::ScalarImageToListAdaptor< ImageType > AdaptorType;
			AdaptorType::Pointer adaptor = AdaptorType::New();
			adaptor->SetImage( reader->GetOutput() );

			typedef PixelType HistogramMeasurementType;
			typedef itk::Statistics::ListSampleToHistogramGenerator <AdaptorType, HistogramMeasurementType > GeneratorType;
			GeneratorType::Pointer generator = GeneratorType::New();
			
			typedef GeneratorType::HistogramType HistogramType;
			HistogramType::SizeType size;
			int nOfBins = 255;
			size.Fill(nOfBins);
			generator->SetListSample(adaptor);
			generator->SetNumberOfBins(size);
			generator->SetMarginalScale(10.0);
			generator->SetHistogramMin(calculator->GetMinimum());
			generator->SetHistogramMax(calculator->GetMaximum());
			generator->Update();
			
			HistogramType::ConstPointer histogram = generator->GetOutput();
			const unsigned int histogramSize = histogram->Size();
			int interval = (calculator->GetMaximum() - calculator->GetMinimum()) / (nOfBins-1);
			
			//std::cout << "Histogram: " << histogramSize << std::endl;
			//std::cout << "x	y"  << std::endl;
			
			// entropy calculation
			// without background
			int countOfNonZeros = 0;
			for( unsigned int bin=0; bin < histogramSize; bin++ ){
				//std::cout << calculator->GetMinimum() + (signed)bin*interval << ' ' << histogram->GetFrequency( bin, 0 ) << std::endl;	
				if (bin != 0 && histogram->GetFrequency( bin, 0 ) != 0) // excluding the first bin (background) and the zero frequencies
					countOfNonZeros ++;
			}
			
			vnl_matrix<double> histo;
			histo.set_size(countOfNonZeros,2);
			countOfNonZeros = 0;
			for( unsigned int bin=0; bin < histogramSize; bin++ ){
				if (bin != 0 && histogram->GetFrequency( bin, 0 ) != 0){ // excluding the first bin (background) and the zero frequencies
					histo(countOfNonZeros,0) = calculator->GetMinimum() + (signed)bin*interval;
					histo(countOfNonZeros,1) = histogram->GetFrequency( bin, 0 );
					countOfNonZeros ++;
				}
			}
			/*
			// with background
			vnl_matrix<double> histo;
			histo.set_size(nOfBins,2);
			for( unsigned int bin=0; bin<nOfBins; bin++ ){
				histo(bin,0) = calculator->GetMinimum() + (signed)bin*interval;
				histo(bin,1) = histogram->GetFrequency( bin, 0 );
			}
*/
			std::cout << "writing histogram" << std::endl;
			VnlWriterMatrix* writer = new VnlWriterMatrix;
			QString temp = _imageFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
			}
			temp.append("_histogram.txt");
			writer->SetFileName(temp);
			writer->SetVnlMatrix(histo);
			writer->MatrixShapeUpdate();

			std::cout << "calculating the entropy" << std::endl;
			StatisticsEntropy* statisticsEntropy = new StatisticsEntropy;
			vnl_vector<double> histoFreq;
			histoFreq.set_size(histo.rows());
			histoFreq = histo.get_column(1);
			statisticsEntropy->SetVector(histoFreq);
			statisticsEntropy->CalculateShannonEntropy();
			std::cout << statisticsEntropy->GetEntropy() << std::endl;
			entropy(i) = statisticsEntropy->GetEntropy();

			
			// cleaning
			delete writer;
			delete statisticsEntropy;
		}

		VnlWriterVector* writer = new VnlWriterVector;
		QString temp = _imageFileNames[0];
		if (temp.lastIndexOf("/") == -1){
			temp.remove(temp.lastIndexOf("\\")+1,temp.size()-1);
		}
		else {
			temp.remove(temp.lastIndexOf("/")+1,temp.size()-1);
		}
		temp.append("entropy.txt");
		writer->SetFileName(temp);
		writer->SetVnlVector(entropy);
		writer->Update();

		std::cout << "entropy file written" << std::endl;

	}

	void ImageHandler::CombineMasks(){

		for (int i=0; i<_maskFileNamesTwo.size(); i++){
		
			// load the three-layer mask
			ImageHandler* reader = new ImageHandler;
			std::cout << "three-layer mask: " << _maskFileNamesTwo[i].ascii() << std::endl;
			reader->SetImageFileName(_maskFileNamesTwo[i].ascii());
			reader->MetafileReaderUpdate();

			// load the correspondent binary mask
			ImageHandler* readerTwo = new ImageHandler;
			QString temp = _maskFileNamesTwo[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("\\")+1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("/")+1);
			}
			for (int y=0; y<_maskFileNames.size(); y++){
				if (_maskFileNames[y].contains(temp)){
					std::cout << "binary mask: " << _maskFileNames[y].ascii() << std::endl;
					readerTwo->SetImageFileName(_maskFileNames[y]);
					readerTwo->MetafileReaderUpdate();
				}
			}
			
			// load the correspondent segmented image to get the right image info 
			// (in some cases the segmented and classified images have different sizes) 
			
			// new combined mask
			ImageHandler::ImageType::Pointer newMask = ImageType::New();
						
			temp = _maskFileNamesTwo[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("\\")+1);
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.remove(0, temp.lastIndexOf("/")+1);
			}
			for (int y=0; y<_imageFileNames.size(); y++){
				if (_imageFileNames[y].contains(temp)){
									
					// uncomprehensible commands from itk
					typedef itk::ImageSeriesReader< ImageHandler::ImageType > ReaderType;
					ReaderType::Pointer reader = ReaderType::New();
					
					typedef itk::GDCMImageIO ImageIOType;
					ImageIOType::Pointer dicomIO = ImageIOType::New();
					reader->SetImageIO( dicomIO );

					typedef itk::GDCMSeriesFileNames NamesGeneratorType;
					NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
					nameGenerator->SetUseSeriesDetails( true );
					std::cout << "segmented: " << _imageFileNames[y].ascii() << std::endl;
					nameGenerator->SetDirectory( _imageFileNames[y].ascii() );
					
					typedef std::vector< std::string > SeriesIdContainer;
					const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
					SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
					SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
					while( seriesItr != seriesEnd )
					{
					//std::cout << seriesItr->c_str() << std::endl;
					seriesItr++;
					}
					std::string seriesIdentifier;
					seriesIdentifier = seriesUID.begin()->c_str();

					typedef std::vector< std::string > FileNamesContainer;
					FileNamesContainer fileNames;
					fileNames = nameGenerator->GetFileNames( seriesIdentifier );
					typedef itk::ImageSeriesReader< ImageHandler::ImageType > ImageSeriesReaderType;
					ImageSeriesReaderType::Pointer imageSeriesReader = ImageSeriesReaderType::New();
					imageSeriesReader->SetFileNames( fileNames );
					imageSeriesReader->Update();

					// image data
					const ImageHandler::ImageType::SpacingType& spacing = imageSeriesReader->GetOutput()->GetSpacing();
					const ImageHandler::ImageType::PointType& origin = imageSeriesReader->GetOutput()->GetOrigin();
					const ImageHandler::ImageType::SizeType& size = imageSeriesReader->GetOutput()->GetBufferedRegion().GetSize();
					std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
					std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
					std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
					
					newMask->SetRegions(imageSeriesReader->GetOutput()->GetRequestedRegion() );
					newMask->CopyInformation( imageSeriesReader->GetOutput() );
					
				}
			}

			// new combined mask
			newMask->Allocate();
			newMask->FillBuffer(0);
			const ImageHandler::ImageType::SpacingType& spacing = newMask->GetSpacing();
			const ImageHandler::ImageType::PointType& origin = newMask->GetOrigin();
			const ImageHandler::ImageType::SizeType& size = newMask->GetBufferedRegion().GetSize();
			std::cout << "new image" << std::endl;
			std::cout << "spacing: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << std::endl;
			std::cout << "origin: " << origin[0] << ' ' << origin[1] << ' ' << origin[2] << std::endl;
			std::cout << "size: " << size[0] << ' ' << size[1] << ' ' << size[2] << std::endl;
					

			// iterators
			typedef itk::ImageRegionIterator< ImageType> IteratorType;
			IteratorType it1 (reader->GetImage(), reader->GetImage()->GetRequestedRegion());
			IteratorType it2 (readerTwo->GetImage(), readerTwo->GetImage()->GetRequestedRegion());
			IteratorType it3 (newMask, newMask->GetRequestedRegion());

			// create image
			for (it2.GoToBegin(),it3.GoToBegin(); !it2.IsAtEnd(); ++it2,++it3){
				if (it2.Get() == 1)
					it3.Set(300); // cortical
			}
			for (it1.GoToBegin(),it3.GoToBegin(); !it1.IsAtEnd(); ++it1,++it3){
				if (it1.Get() == 200)
					it3.Set(200); // trabecular
				else if(it1.Get() == 100)
					it3.Set(100); // marrow
			}
		
			// write image
			ImageHandler* writer = new ImageHandler;
			QString combinedMaskFileName = ("C:/0.Data/0. Original data/3. femur left masks/3. combination/");
			QString tempTwo = _maskFileNamesTwo[i];
			if (tempTwo.lastIndexOf("/") == -1){
				tempTwo.remove(0, tempTwo.lastIndexOf("\\")+1);
			}
			else {
				tempTwo.remove(0, tempTwo.lastIndexOf("/")+1);
			}
			combinedMaskFileName.append(tempTwo);
			std::cout << "combinedMask: " << combinedMaskFileName.ascii() << std::endl;
			writer->SetImageFileName(combinedMaskFileName);
			writer->SetImage(newMask);
			writer->MetafileWriterUpdate();

			delete reader;
			delete readerTwo;
			delete writer;

		}
	}

	void ImageHandler::ExtractMasks(){
		
		// read mask to split
		_maskFileName = "C:/0.Data/test data/extract masks - ben/s339.mhd";
		ImageHandler* reader = new ImageHandler;
		std::cout << "mask: " << _maskFileName.toAscii().data() << std::endl;
		reader->SetImageFileName(_maskFileName);
		reader->MetafileReaderUpdate();

		// new mask
		ImageHandler::ImageType::Pointer newMask = ImageType::New();
		newMask->SetRegions(reader->GetImage()->GetRequestedRegion() );
		newMask->CopyInformation( reader->GetImage());
		newMask->Allocate();
		newMask->FillBuffer(0);

		// iterators
		typedef itk::ImageRegionIterator< ImageType> IteratorType;
		IteratorType it1 (reader->GetImage(), reader->GetImage()->GetRequestedRegion());
		IteratorType it2 (newMask, newMask->GetRequestedRegion());
		
		// create image
		for (it1.GoToBegin(),it2.GoToBegin(); !it1.IsAtEnd(); ++it1,++it2){
			if (it1.Get() == 2)
				it2.Set(1); 
		}
		
		// write image
		ImageHandler* writer = new ImageHandler;
		writer->SetImageFileName("cement_s339.mhd");
		writer->SetImage(newMask);
		writer->MetafileWriterUpdate();





		
	
	}
	void ImageHandler::CreateTestImages(){
	
		/*
		// create scalar image
		typedef signed short VoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Image< VoxelType, Dimension > ImageType;
		ImageType::Pointer image = ImageType::New();
		ImageType::IndexType start;
		start[0] = 0; // first index on X
		start[1] = 0; // first index on Y
		start[2] = 0; // first index on Z
		ImageType::SizeType size;
		size[0] = 2; // size along X
		size[1] = 2; // size along Y
		size[2] = 2; // size along Z
		ImageType::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );
		image->SetRegions( region );
		image->Allocate();

		ImageType::IndexType pixelIndex;
		ImageType::PixelType pixelValue;

		pixelIndex[0] = 0; pixelIndex[1] = 0; pixelIndex[2] = 0; pixelValue = 301; image->SetPixel( pixelIndex, pixelValue);
		pixelIndex[0] = 1; pixelIndex[1] = 0; pixelIndex[2] = 0; pixelValue = 302; image->SetPixel( pixelIndex, pixelValue);
		pixelIndex[0] = 0; pixelIndex[1] = 1; pixelIndex[2] = 0; pixelValue = 303; image->SetPixel( pixelIndex, pixelValue);
		pixelIndex[0] = 1; pixelIndex[1] = 1; pixelIndex[2] = 0; pixelValue = 304; image->SetPixel( pixelIndex, pixelValue);
		pixelIndex[0] = 0; pixelIndex[1] = 0; pixelIndex[2] = 1; pixelValue = 305; image->SetPixel( pixelIndex, pixelValue);
		pixelIndex[0] = 1; pixelIndex[1] = 0; pixelIndex[2] = 1; pixelValue = 306; image->SetPixel( pixelIndex, pixelValue);
		pixelIndex[0] = 0; pixelIndex[1] = 1; pixelIndex[2] = 1; pixelValue = 307; image->SetPixel( pixelIndex, pixelValue);
		pixelIndex[0] = 1; pixelIndex[1] = 1; pixelIndex[2] = 1; pixelValue = 308; image->SetPixel( pixelIndex, pixelValue);

		typedef itk::ImageRegionIterator< ImageType> IteratorType;	
		IteratorType it1 (image, image->GetRequestedRegion());
		
		for (it1.GoToBegin(); !it1.IsAtEnd(); ++it1)
			std::cout << it1.Get() << std::endl;

		typedef itk::ImageFileWriter< ImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName("Intensity1.mhd");
		writer->SetInput( image );
		writer->Update();
		std::cout << "metafile written" << std::endl;
		*/

		// create vector image
		typedef float FieldVoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Vector<FieldVoxelType,Dimension> VectorType;
		typedef itk::Image<VectorType,Dimension> FieldType; 
		FieldType::Pointer field = FieldType::New();
		FieldType::IndexType start;
		start[0] = 0; // first index on X
		start[1] = 0; // first index on Y
		start[2] = 0; // first index on Z
		FieldType::SizeType size;
		size[0] = 2; // size along X
		size[1] = 2; // size along Y
		size[2] = 2; // size along Z
		FieldType::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );
		field->SetRegions( region );
		field->Allocate();

		FieldType::IndexType pixelIndex;
		FieldType::PixelType pixelValue;

		pixelIndex[0] = 0; pixelIndex[1] = 0; pixelIndex[2] = 0; 
		pixelValue[0] = 80; pixelValue[1] = 82; pixelValue[2] = 84; 
		field->SetPixel( pixelIndex, pixelValue);
		
		pixelIndex[0] = 1; pixelIndex[1] = 0; pixelIndex[2] = 0; 
		pixelValue[0] = 86; pixelValue[1] = 88; pixelValue[2] = 90; 
		field->SetPixel( pixelIndex, pixelValue);
		
		pixelIndex[0] = 0; pixelIndex[1] = 1; pixelIndex[2] = 0; 
		pixelValue[0] = 92; pixelValue[1] =94; pixelValue[2] = 96; 
		field->SetPixel( pixelIndex, pixelValue);
		
		pixelIndex[0] = 1; pixelIndex[1] = 1; pixelIndex[2] = 0; 
		pixelValue[0] = 98; pixelValue[1] = 100; pixelValue[2] = 102; 
		field->SetPixel( pixelIndex, pixelValue);
		
		pixelIndex[0] = 0; pixelIndex[1] = 0; pixelIndex[2] = 1; 
		pixelValue[0] = 104; pixelValue[1] = 106; pixelValue[2] = 108; 
		field->SetPixel( pixelIndex, pixelValue);
		
		pixelIndex[0] = 1; pixelIndex[1] = 0; pixelIndex[2] = 1; 
		pixelValue[0] = 110; pixelValue[1] = 112; pixelValue[2] = 114; 
		field->SetPixel( pixelIndex, pixelValue);
		
		pixelIndex[0] = 0; pixelIndex[1] = 1; pixelIndex[2] = 1; 
		pixelValue[0] = 116; pixelValue[1] = 118; pixelValue[2] = 120; 
		field->SetPixel( pixelIndex, pixelValue);
		
		pixelIndex[0] = 1; pixelIndex[1] = 1; pixelIndex[2] = 1; 
		pixelValue[0] = 122; pixelValue[1] = 124; pixelValue[2] =126; 
		field->SetPixel( pixelIndex, pixelValue);

				
		typedef itk::ImageRegionIterator< FieldType> IteratorType;	
		IteratorType it1 (field, field->GetRequestedRegion());
		
		for (it1.GoToBegin(); !it1.IsAtEnd(); ++it1){
			std::cout << it1.Get()[0] << std::endl;
			std::cout << it1.Get()[1] << std::endl;
			std::cout << it1.Get()[2] << std::endl;
		}

		typedef itk::ImageFileWriter< FieldType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName("00003.mhd");
		writer->SetInput( field );
		writer->Update();
		std::cout << "metafile written" << std::endl;



	

		


		
	}

	void ImageHandler::SubtractImages(){

		// new image
		ImageType::Pointer outputImage = ImageType::New();
		outputImage->SetRegions(_image->GetRequestedRegion() );
		outputImage->CopyInformation( _image);
		outputImage->Allocate();
		outputImage->FillBuffer(0);
		
		// subtraction
		typedef itk::ImageRegionIterator< ImageType> IteratorType;	
		IteratorType inputIt( _image, _image->GetRequestedRegion() );
		IteratorType inputIt2( _image2, _image2->GetRequestedRegion() );
		IteratorType outputIt( outputImage, outputImage->GetRequestedRegion() );
		
		for ( inputIt.GoToBegin(), inputIt2.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++inputIt2, ++outputIt){
			
			outputIt.Set( inputIt.Get() - inputIt2.Get());
		}

		// saving
		typedef itk::ImageFileWriter< ImageType > ImageWriterType;
		ImageWriterType::Pointer imageWriter = ImageWriterType::New();
		imageWriter->SetFileName("subtracted.mhd");	
		imageWriter->SetInput(outputImage);
		imageWriter->Update();
		

	
	}
}


/*void ImageHandler::IntensitiesToVector{

		for (int i=0; i<_imageFileNames.size(); i++){
			
			// reading the image
			typedef itk::ImageFileReader< ImageHandler::ImageType > ReaderType;
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(_imageFileNames[i].ascii());
			std::cout <<_imageFileNames[i].ascii() << std::endl;
			try{
				reader->Update();
			}
			catch( itk::ExceptionObject & excp ){
				std::cerr << "!!! Problem reading the input file" << std::endl;
			}
			
			typedef itk::MinimumMaximumImageCalculator< ImageHandler::ImageType > MinMaxCalculator;
			MinMaxCalculator::Pointer calculator = MinMaxCalculator::New();
			calculator->SetImage(reader->GetOutput());
			calculator->Compute();
			std::cout << "min intensity: " << calculator->GetMinimum() << "; max intensity: " << calculator->GetMaximum() << std::endl; 
		
			// calculating the number of non-background voxels
			ImageType::IndexType pixelIndex;
			pixelIndex[0] = 2; pixelIndex[1] = 2; pixelIndex[2] = 2;
			ImageType::PixelType backgroundValue = reader->GetOutput()->GetPixel( pixelIndex );
			typedef itk::ImageRegionIterator< ImageType> IteratorType;
			IteratorType it( reader->GetOutput(), reader->GetOutput()->GetRequestedRegion() );
			int count = 0;
			for ( it.GoToBegin(); !it.IsAtEnd(); ++it){
				if (it.Get() != backgroundValue)
					count ++;
			}

			// putting the non-background voxels

		
		}
	

	}
	*/
	


	