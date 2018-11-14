/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshMorpherPointPropagation.h>

#include <VectorImageHandler.h>
#include <PointReaderXyzId.h>
#include <PointWriterXyz.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>//

#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkSTLReader.h>

#include <QString>
#include <QFile>
#include <QString>
#include <QTextStream>


using namespace points;
using namespace image;
using namespace mesh;

namespace mesh{

	// constructor
	MeshMorpherPointPropagation::MeshMorpherPointPropagation(){
	}

	// destructor
	MeshMorpherPointPropagation::~MeshMorpherPointPropagation(){
	}

	// member functions
	void MeshMorpherPointPropagation::Update(){
			
		// filenames
		QString txtFileName;
		for (int i=0; i<_fileNames.size(); i++){
			if (_fileNames[i].endsWith(".txt")){
				txtFileName=_fileNames[i].ascii();
				_fileNames.removeAt(i);
			}
		}
	
		// reference: read .txt file and points in vtkPoints
		vtkPoints* referenceLandmarks = vtkPoints::New();
		vtkDoubleArray* referenceID = vtkDoubleArray::New();
		PointReaderXyzId* pointIdReader = new PointReaderXyzId;
		pointIdReader->SetFileName(txtFileName);
		pointIdReader->Update();
		referenceLandmarks = pointIdReader->GetOutput();
		referenceID = pointIdReader->GetIdVector();
				
		// for each VF 
		for (int i=0; i<_fileNames.size(); i++){
			
			// load VF
			std::cout << "load VF" <<  std::endl;
			std::cout << _fileNames[i].ascii() << std::endl;
			const unsigned int VectorDimension = 3;
			typedef itk::Vector< float, VectorDimension > PixelType;
			const unsigned int ImageDimension = 3;
			typedef itk::Image< PixelType, ImageDimension > FieldType;
			typedef itk::ImageFileReader< FieldType > FieldReaderType;
			FieldReaderType::Pointer reader = FieldReaderType::New();
			reader->SetFileName(_fileNames[i].ascii());
			try
			{
				reader->Update();
			}
			catch( itk::ExceptionObject & excp )
			{
				std::cerr << excp << std::endl;
			}

			// VF to DVF
			std::cout << "convert to DVF" <<  std::endl;
			VectorImageHandler* imageHandler = new VectorImageHandler;
			imageHandler->SetField(reader->GetOutput());
			imageHandler->VFtoDVFinverted();
			imageHandler->ITKtoVTKtoITK();

			// extract correspondent moving points 
			vtkPoints* movingPoints = vtkPoints::New();
			const FieldType::SpacingType& spacing = imageHandler->GetField()->GetSpacing();
			 
			for (int a=0; a<referenceLandmarks->GetNumberOfPoints(); a++){
				double refCoord[3]; 
				double movCoord[3];
				
				referenceLandmarks->GetPoint(a, refCoord);
				FieldType::IndexType pixelIndex;
				
				pixelIndex[0] = std::floor(refCoord[0] / spacing[0] + 0.5);
				pixelIndex[1] = std::floor(refCoord[1] / spacing[1] + 0.5);
				pixelIndex[2] = std::floor(refCoord[2] / spacing[2] + 0.5);
				
				FieldType::PixelType pixelValue = imageHandler->GetField()->GetPixel( pixelIndex );
				
				movCoord[0] =  refCoord[0] + pixelValue[0];
				movCoord[1] =  refCoord[1] + pixelValue[1];
				movCoord[2] =  refCoord[2] + pixelValue[2];
							
				movingPoints->InsertPoint(a, movCoord); 
			}
			
			// load .stl
			_fileNames[i].replace(QString("_velocity_field.mhd"), QString("_affine.stl"));
			std::cout << "load mesh " << std::endl;
			std::cout << _fileNames[i].ascii() << std::endl;
			vtkSTLReader* stlReader = vtkSTLReader::New();
			stlReader->SetFileName(_fileNames[i].ascii());
			stlReader->Update();

			// find moving landmarks from the .stl
			vtkPoints* movingLandmarks = vtkPoints::New();
			for (int a=0; a<movingPoints->GetNumberOfPoints(); a++){
				
				double coord[3];
				movingPoints->GetPoint(a, coord);
				int index = 0;
				findClosestNode (coord, stlReader->GetOutput()->GetPoints(), index);
				stlReader->GetOutput()->GetPoint(index,coord);
				movingLandmarks->InsertPoint(a, coord);
			
			}
		
			// write .txt file for moving landmarks (to check their position visually)
			_fileNames[i].replace(QString(".stl"), QString(".txt"));
			PointWriterXyz* pointWriter = new PointWriterXyz;
			pointWriter->SetFileName(_fileNames[i].ascii());
			pointWriter->SetInput(movingLandmarks);
			pointWriter->Update();
		
			// write .xml file
			_fileNames[i].replace(QString(".txt"), QString(".xml"));
			writeXmlFile(_fileNames[i].ascii(), movingLandmarks, referenceID);
			std::cout << "write .xml file" <<  std::endl;
		}
	
	}


	
	// private functions
	void MeshMorpherPointPropagation::findClosestNode (double point[3], vtkPoints* points, int &minIndex){
	
		// find the absolute distances between the given application point and all the mesh nodes
		// inputs:
		// - point to loot for the closest to
		// - point cloud (mesh nodes)
		// output:
		// - closest point id

		double min = 0.0; //int minIndex;
		int nOfCoord = 3;
		
		for (int i=0; i<points->GetNumberOfPoints(); i++){
			double coord[3];
			points->GetPoint(i, coord);

			double distance = 0.0;
			for (int a=0; a<3; a++)
				distance += ((coord[a]-(point[a])) * (coord[a]-(point[a])));
				distance = sqrt(distance);
			if (i==0){
				min = distance;
				minIndex=i+1; // in abaqus node id starts from 1
			}
			else{
				if (distance < min){
					min = distance;
					minIndex=i+1; // in abaqus node id starts from 1
					}
				}
		}
	
	}
	void MeshMorpherPointPropagation::writeXmlFile(QString fileName, vtkPoints* movingLandmarks, vtkDoubleArray* referenceID){
	
		// ! change the reference file name here 
		QString refCdbFileName = ("00089.cdb");

		QFile outFileXml(fileName.toAscii().data());
		outFileXml.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFileXml(&outFileXml);
		fileName.replace(QString(".xml"), QString(".stl"));
		fileName.remove(0,fileName.lastIndexOf("/")+1);
		double coord[3];
		writeFileXml << "<morphing>" << endl;
		writeFileXml << "	<template_mesh type=\"vol\">.\\\\..\\\\test\\\\" << refCdbFileName << "</template_mesh>" << endl;
		writeFileXml << "	<STL_file coef=\"1.\">.\\\\..\\\\test\\\\" << fileName.toAscii().data() << "</STL_file>" << endl;
		writeFileXml << "   <same_femur_orientation>true</same_femur_orientation>" << endl;
		writeFileXml << "   <debug_level>2</debug_level>" << endl;
		writeFileXml << "   <anatomical_points_stl method =\"coord\">" << endl; // moving landmarks
		movingLandmarks->GetPoint(0, coord);
		writeFileXml << "		<AP id=\"1\" id_name=\"head\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(1, coord);		
		writeFileXml << "		<AP id=\"2\" id_name=\"big_trochanter\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(2, coord);
		writeFileXml << "		<AP id=\"3\" id_name=\"under_big_trochanter\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1]<< "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(3, coord);
		writeFileXml << "		<AP id=\"4\" id_name=\"small_trochanter\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(4, coord);
		writeFileXml << "		<AP id=\"5\" id_name=\"medial_condyle\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(5, coord);
		writeFileXml << "		<AP id=\"6\" id_name=\"lateral_condyle\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(6, coord);
		writeFileXml << "		<AP id=\"7\" id_name=\"down_medial_side\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(7, coord);
		writeFileXml << "		<AP id=\"8\" id_name=\"down_lateral_side\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(8, coord);
		writeFileXml << "		<AP id=\"9\" id_name=\"down_behind\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		movingLandmarks->GetPoint(9, coord);
		writeFileXml << "		<AP id=\"10\" id_name=\"middle\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		writeFileXml << "	</anatomical_points_stl>" << endl;
		writeFileXml << "	<anatomical_points_genmesh method =\"node\">" << endl;
		int id;
		id=referenceID->GetValue(0);
		writeFileXml << "		<AP id=\"1\" id_name=\"head\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(1);
		writeFileXml << "		<AP id=\"2\" id_name=\"big_trochanter\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(2);
		writeFileXml << "		<AP id=\"3\" id_name=\"under_big_trochanter\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(3);
		writeFileXml << "		<AP id=\"4\" id_name=\"small_trochanter\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(4);
		writeFileXml << "		<AP id=\"5\" id_name=\"medial_condyle\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(5);
		writeFileXml << "		<AP id=\"6\" id_name=\"lateral_condyle\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(6);
		writeFileXml << "		<AP id=\"7\" id_name=\"down_medial_side\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(7);
		writeFileXml << "		<AP id=\"8\" id_name=\"down_lateral_side\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(8);
		writeFileXml << "		<AP id=\"9\" id_name=\"down_behind\"> " << id << " </AP>" << endl;	
		id=referenceID->GetValue(9);
		writeFileXml << "		<AP id=\"10\" id_name=\"middle\"> " << id << " </AP>" << endl;	
		writeFileXml << "	</anatomical_points_genmesh>" << endl;
		writeFileXml << "	<resdir>.\\\\..\\\\test\\\\results</resdir>" << endl;
		writeFileXml << "	<res_surface_mesh>SurfaceResult.cdb</res_surface_mesh>" << endl;
		writeFileXml << "	<res_volume_mesh>VolumeResult.cdb</res_volume_mesh>" << endl;
		writeFileXml << "</morphing>" << endl;
		outFileXml.close();
		
	}

}

