/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef IMAGEHANDLER_H
#define IMAGEHANDLER_H

#include <itkImage.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <QString>
#include <QStringList>


namespace image{

	/** 
	* It executes many operations on images that are commonly useful
	*
	* ImageType is defined as unsigned short of dimension 3.
	*/
	class ImageHandler{
		
	public:
		
		// 2D image
		typedef float PixelType; // used in the image combined pca
		static const unsigned int PixelDimension = 2;
		typedef itk::Image< PixelType, PixelDimension > PixelImageType;

		// 3D image
		typedef signed short VoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Image< VoxelType, Dimension > ImageType;

		// vector image
		typedef float FieldVoxelType;
		typedef itk::Vector<FieldVoxelType,Dimension> VectorType;
		typedef itk::Image<VectorType,Dimension> FieldType; 


		// constructor/destructor
		ImageHandler();
		~ImageHandler();

		// accessors
		/**
		* Sets the dicom series folder name
		*/
		void SetDicomSerieFolderName (QString imageFolderFileName) {_imageFolderFileName = imageFolderFileName;}
		/**
		* Sets the dicom series folder name
		*/
		void SetDicomMaskFolderName (QString maskFolderFileName) {_maskFolderFileName = maskFolderFileName;}
		/**
		* Sets the image file name
		*/
		void SetImageFileName (QString imageFileName) {_imageFileName = imageFileName;}
		/**
		* Sets the image file names
		*/
		void SetImageFileNames (QStringList imageFileNames) {_imageFileNames = imageFileNames;}
		/**
		* Sets the vector image file name
		*/
		void SetFieldFileName (QString fieldFileName) {_fieldFileName = fieldFileName;}
		/**
		* Sets the vector image file names
		*/
		void SetFieldFileNames (QStringList fieldFileNames) {_fieldFileNames = fieldFileNames;}
		/**
		* Sets the mask file name
		*/
		void SetMaskFileName (QString maskFileName) {_maskFileName = maskFileName;}
		/**
		* Sets the mask file names
		*/
		void SetMaskFileNames (QStringList maskFileNames) {_maskFileNames = maskFileNames;}
		/**
		* Sets the mask file names 2
		*/
		void SetMaskFileNamesTwo (QStringList maskFileNamesTwo) {_maskFileNamesTwo = maskFileNamesTwo;}
		/**
		* Sets the image
		*/
		void SetImage (ImageType::Pointer image) {_image = image;}
		/**
		* Sets the image 2
		*/
		void SetImage2 (ImageType::Pointer image) {_image2 = image;}
		/**
		* Sets the threshold for the marching cubes. Default threshold = -1000.
		*/
		void SetThreshold (int threshold) {_threshold = threshold;}
		/**
		* Sets the voxel size. Default spacing = 1.25 x 1.25 x 1.25.
		*/
		void SetSpacing (ImageType::SpacingType spacing) {_spacing = spacing;}
		/**
		* Sets the vector image
		*/
		void SetDeformationVectorField(FieldType::ConstPointer constField) {_constField = constField;}
		/**
		* Gets the image
		*/
		ImageType::Pointer GetImage () {return _image;}
		/**
		* Gets the field
		*/
		FieldType::Pointer GetField () {return _field;}
		/**
		* Gets the mesh
		*/
		vtkPolyData* GetMesh() {return _mesh;}
		/**
		* Gets the vtk image
		*/
		vtkImageData* GetVtkImage() {return _vtkImage;}

		
		// member functions
		
		/******************************** IMAGE READERS ********************************/
		/**
		* Reads dicom series. Use SetDicomSerieFolderName() and GetImage()
		*/
		void DicomReaderUpdate();
		/**
		* Reads dicom series. Use SetDicomMaskFolderName() and GetImage()
		*/
		void DicomMaskReaderUpdate();
		/**
		* Reads metafiles. Use SetImageFileName() and GetImage()
		*/
		void MetafileReaderUpdate();
		/**
		* Reads mask images. Use SetMaskFileName()
		*/
		void MaskReaderUpdate();
		/**
		* Reads vector field. Use SetFieldFileName()
		*/
		void FieldReaderUpdate();

		
		
		/******************************** IMAGE WRITERS ********************************/
		/**
		* Writes metafiles. Use SetImageFileName() and SetImage()
		*/
		void MetafileWriterUpdate();
		/**
		* Writes dicom series. Use SetDicomSerieFolderName() and SetImage(). ! Negative value are not written (possible problems in the itk class)
		*/
		void DicomWriterUpdate();
		
		
		/******************************** IMAGE MASKING ********************************/
		/**
		* Masks an image overwriting the original one. Use SetImageFileName() and SetMaskFileName()
		*/
		void MaskToImage();
		/**
		* Sets the mask grey levels to 0, 100, 200 and 300, without greylevels in between. Use SetImage()
		*/
		void MaskGreyLevelsUpdate();

		
		/******************************** ITK VTK CONVERSION ***************************/
		/**
		* Tranforms an ITK image into a VTK image. Use SetImage(), MetaImageReader() and GetVtkImage()
		*/
		void ITKtoVTK();
		/**
		* Tranforms an ITK image into a VTK image and then back to an ITK image.
		* This is done to have images and meshes in the same coordinate systems
		* Use SetImage(), MetaImageReader()
		*/
		void ITKtoVTKtoITK();
		
		
		/*********************************** PROCESSING ********************************/
		/**
		* Transforms an image in a mesh. Use SetImage(), MetaImageReader(), SetThreshold(), and GetMesh()
		*/
		void MarchingCubesUpdate();
		/**
		* Resemples the image. Use SetImage(), SetSpacing and GetImage()
		*/
		void ResampleUpdate();
		/**
		* Applies the DVF to an image and warps it. Use SetDVF(), SetImage() and GetImage()
		*/
		void ApplyDVFToImageUpdate();
		/**
		* Averages grey level of images. Use SetImageFileNames()
		*/
		void AverageImageGreyLevels();
		/**
		* Executes a fake calibration on CT images. 0 HU are put to 0.9937mg/mm3 (water), the max HU is put to 1.2mg/mm3 (compact bone)
		* !!! All values are multiplied by 1000 in order to have short images (could be considered as g/mm3)
		*/
		// void PseudoCalibration(); (moved to ImageHandler.h)
		/**
		* Extrudes the image. First the voxels close to the background are deleted in order to avoid the partial volume effect.
		* Then the new outer voxels are averaged and put outside for 3 times
		*/
		void Extrusion();
		/**
		* Crops and existing image, creating a black frame around it. Dimensions are hardcoded.
		*/
		void Crop();
		/**
		* Reads the dicom field that are currently set from the first slice of each stack.
		*/
		void DicomHeaderReaderFromSlice();
		/**
		* Calculates the image histogram.
		*/
		void ImageHistogram();
		/**
		* Combines the mask images (my segmentation with stryker one). 
		* Use SetMaskFileNames() for the binary ones and SetMaskFileNamesTwo() for the stryker ones. 
		* Set manually in the function the output folder
		*/
		void CombineMasks();
		/**
		* Extracts the layers of the mask and write it to different images
		*/
		void ExtractMasks();
		/**
		* Creates test images for debugging
		*/
		void CreateTestImages();
		/**
		* Subtract _image and _image2
		*/
		void SubtractImages();
					

	protected:
		
		// data members from accessors
		QString _imageFolderFileName;
		QString _maskFolderFileName;
		QString _imageFileName;
		QStringList _imageFileNames;
		QString _fieldFileName;
		QStringList _fieldFileNames;
		QString _maskFileName;
		QStringList _maskFileNames;
		QStringList _maskFileNamesTwo;
		ImageType::Pointer _image;
		ImageType::Pointer _image2;
		ImageType::Pointer _mask;
		vtkPolyData* _mesh;
		vtkImageData* _vtkImage;
		int _threshold;
		ImageType::SpacingType _spacing;
		FieldType::Pointer _field;
		FieldType::ConstPointer _constField;
		int _edgePaddingValue;
		
	};
}
#endif //IMAGEHANDLER_H