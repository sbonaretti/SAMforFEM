/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef IMAGEHANDLERFLOAT_H
#define IMAGEHANDLERFLOAT_H

#include <itkImage.h>
#include <itkCenteredAffineTransform.h>
#include <vtkImageData.h>
#include <vtkPoints.h>
#include <QString>

#include <VnlWriterVector.h>

namespace image{

	/** 
	* ImageType is defined as float of dimension 3.
	*/
	class ImageHandlerFloat{
		
	public:
		
		// 3D image
		typedef float VoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Image< VoxelType, Dimension > ImageType;

		// vector image
		typedef float FieldVoxelType;
		typedef itk::Vector<FieldVoxelType,Dimension> VectorType;
		typedef itk::Image<VectorType,Dimension> FieldType; 

		typedef itk::CenteredAffineTransform< double,3 > AffineTransformType;
		
		// constructor/destructor
		ImageHandlerFloat();
		~ImageHandlerFloat();

		// accessors
		/**
		* Sets the image file name
		*/
		void SetImageFileName (QString imageFileName) {_imageFileName = imageFileName;}
		/**
		* Sets the image
		*/
		void SetImage (ImageType::Pointer image) {_image = image;}
		/**
		* Sets the mask
		*/
		void SetMask (ImageType::Pointer mask) {_mask = mask;}
		/**
		* Sets the velocityField (constPointer in order to be compatible with VectorImageHandler.h)
		*/
		void SetSVF (FieldType::ConstPointer svf) {_svf = svf;}
		/**
		* Gets the image
		*/
		ImageType::Pointer GetImage () {return _image;}
		/**
		* Sets the points of the rotation plain
		*/
		void SetPoints (vtkPoints* points) { _points = points;}
		/**
		* Gets the vnl vector
		*/
		vnl_vector<double> GetVector () {return _vector;}
		

		
		// member functions
		/**
		* Reads metafiles. Use SetImageFileName() and GetImage()
		*/
		void MetafileReaderUpdate();
		/**
		* Writes metafiles. Use SetImageFileName() and SetImage()
		*/
		void MetafileWriterUpdate();
		/**
		* Find the minimum of the marrow area. Use SetMark() and SetSVF()
		*/
		void FindMinMarrow();
		/**
		* Executes a fake calibration on CT images. 0 HU are put to 0.9937mg/mm3 (water), the max HU is put to 1.2mg/mm3 (compact bone)
		* !!! All values are multiplied by 1000 in order to have short images (could be considered as g/mm3)
		*/
		void PseudoCalibration();
		/**
		* Extrudes the image. First the voxels close to the background are deleted in order to avoid the partial volume effect.
		* Then the new outer voxels are averaged and put outside for 3 times
		*/
		void Extrusion();
		/**
		* Tranforms an ITK image into a VTK image and then back to an ITK image.
		* This is done to have images and meshes in the same coordinate systems
		* Use SetImage(), MetaImageReader()
		*/
		void ITKtoVTKtoITK();
		/**
		* Rotate image to a plain defined by three points
		*/
		void RotateImage();
		/**
		* Extract an area around seed points
		*/
		void ExtractVolumeImage();



	protected:
		
		// data members from accessors
		QString _imageFileName;
		ImageType::Pointer _image;
		ImageType::Pointer _mask;
		FieldType::ConstPointer _svf;
		double _threshold;
		vtkImageData* _vtkImage;
		float _min;
		float _edgePaddingValue;
		vtkPoints* _points;
		vnl_vector<double> _vector;
		
	};

}
#endif //IMAGEHANDLERFLOAT_H
		

