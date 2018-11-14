/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef VECTORIMAGEHANDLER_H
#define VECTORIMAGEHANDLER_H

#include <itkImage.h>

#include <vtkImageData.h>

#include <QString>
#include <QStringList>

namespace image{

	/** 
	* It executes many operations on vector images that are commonly useful
	*
	* ImageType is defined as float of dimension 3.
	*/

	class VectorImageHandler{
		
	public:
		
		// vector image
		static const unsigned int Dimension = 3;
		typedef float FieldVoxelType;
		typedef itk::Vector<FieldVoxelType,Dimension> VectorType;
		typedef itk::Image<VectorType,Dimension> FieldType; 		
		
		// constructor/destructor
		VectorImageHandler();
		~VectorImageHandler();

		// accessors
		/**
		* Sets the vector image file name
		*/
		void SetFieldFileName (QString fieldFileName) {_fieldFileName = fieldFileName;}
		/**
		* Sets the vector image file name list
		*/
		void SetFieldFileNames (QStringList fieldFileNames) {_fieldFileNames = fieldFileNames;}
		/**
		* Sets the vector image
		*/
		void SetField(FieldType::ConstPointer field) {_field = field;}
		/**
		* Gets the vector image
		*/
		FieldType::ConstPointer GetField () {return _field;}
		/**
		* Gets the vtk vector image
		*/
		vtkImageData* GetVtkImageDataOutput() {return _vtkField;}

		// member functions
		
		/**************************** FIELD READER AND WRITER ****************************/
		/**
		* Reads metafiles. Use SetFieldFileName() and GetImage()
		*/
		void MetafileReaderUpdate();
		/**
		* Writes metafiles. Use SetFieldFileName() and SetField()
		*/
		void MetafileWriterUpdate();
		
		/****************************** ITK VTK CONVERSION *******************************/
		/**
		* Tranforms an ITK image into a VTK image.
		*/
		void ITKtoVTK();
		/**
		* Tranforms an ITK image into a VTK image and back to ITK - done to have the same coordinate systems.
		*/
		void ITKtoVTKtoITK();
		
		/******************************* VF DVF CONVERSION *******************************/
		/**
		* Tranforms a velocity field in a deformation vector field to be applied to the reference in order to obtain the moving
		*/
		void VFtoDVFinverted();
		/**
		* Tranforms a velocity field in a deformation vector field to be applied to the moving in order to have the reference shape.
		*/
		void VFtoDVFnotInverted();
		/**
		* Calculates the average velocity field, the absolute values of the distance of each field to the average and prints out the new reference file name.
		* Use SetFieldFileNames()
		*/
		void SelectReferenceAsMinimumDistanceToAverage();
	
	protected:
		// data members from accessors
		QString _fieldFileName;
		QStringList _fieldFileNames;
		FieldType::ConstPointer _field;
		vtkImageData* _vtkField;


	};
}
#endif //VECTORIMAGEHANDLER_H