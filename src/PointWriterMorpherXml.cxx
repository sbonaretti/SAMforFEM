/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <PointWriterMorpherXml.h>

namespace points{

	// constructor
	PointWriterMorpherXml::PointWriterMorpherXml(){
	}
	
	// destructor
	PointWriterMorpherXml::~PointWriterMorpherXml(){
	}

	// overwritten virtual function
	void PointWriterMorpherXml::Update(){

		
		QFile outFileXml(_fileName.toAscii().data());

		outFileXml.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream writeFileXml(&outFileXml);
		_fileName.replace(QString(".xml"), QString(".stl"));
				
		writeFileXml << "<morphing>" << endl;
		writeFileXml << "	<template_mesh type=\"vol\">" << _referenceFileName.toAscii().data() << "</template_mesh>" << endl;
		writeFileXml << "	<STL_file coef=\"1.\">" << _fileName.toAscii().data() << "</STL_file>" << endl;
		writeFileXml << "   <same_femur_orientation>true</same_femur_orientation>" << endl;
		writeFileXml << "   <debug_level>4</debug_level>" << endl;
		writeFileXml << "   <volume_morphing_method>3</volume_morphing_method>" << endl;
		
		writeFileXml << "   <anatomical_points_stl method =\"coord\">" << endl; // moving landmarks
		double coord[3];
		_points->GetPoint(0, coord);
		writeFileXml << "		<AP id=\"1\" id_name=\"head\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(1, coord);		
		writeFileXml << "		<AP id=\"2\" id_name=\"big_trochanter\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(2, coord);
		writeFileXml << "		<AP id=\"3\" id_name=\"under_big_trochanter\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1]<< "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(3, coord);
		writeFileXml << "		<AP id=\"4\" id_name=\"small_trochanter\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(4, coord);
		writeFileXml << "		<AP id=\"5\" id_name=\"medial_condyle\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(5, coord);
		writeFileXml << "		<AP id=\"6\" id_name=\"lateral_condyle\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(6, coord);
		writeFileXml << "		<AP id=\"7\" id_name=\"down_medial_side\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(7, coord);
		writeFileXml << "		<AP id=\"8\" id_name=\"down_lateral_side\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(8, coord);
		writeFileXml << "		<AP id=\"9\" id_name=\"down_behind\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		_points->GetPoint(9, coord);
		writeFileXml << "		<AP id=\"10\" id_name=\"middle\" coordX=\"" << coord[0] << "\" coordY=\"" << coord[1] << "\" coordZ=\"" << coord[2] << "\" ></AP>" << endl;	
		writeFileXml << "	</anatomical_points_stl>" << endl;
		
		writeFileXml << "	<anatomical_points_genmesh method =\"node\">" << endl;
		int id;
		id=_idVector->GetValue(0);
		writeFileXml << "		<AP id=\"1\" id_name=\"head\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(1);
		writeFileXml << "		<AP id=\"2\" id_name=\"big_trochanter\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(2);
		writeFileXml << "		<AP id=\"3\" id_name=\"under_big_trochanter\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(3);
		writeFileXml << "		<AP id=\"4\" id_name=\"small_trochanter\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(4);
		writeFileXml << "		<AP id=\"5\" id_name=\"medial_condyle\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(5);
		writeFileXml << "		<AP id=\"6\" id_name=\"lateral_condyle\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(6);
		writeFileXml << "		<AP id=\"7\" id_name=\"down_medial_side\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(7);
		writeFileXml << "		<AP id=\"8\" id_name=\"down_lateral_side\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(8);
		writeFileXml << "		<AP id=\"9\" id_name=\"down_behind\"> " << id << " </AP>" << endl;	
		id=_idVector->GetValue(9);
		writeFileXml << "		<AP id=\"10\" id_name=\"middle\"> " << id << " </AP>" << endl;	
		writeFileXml << "	</anatomical_points_genmesh>" << endl;
		
		QString temp = _fileName;
	
		if (_fileName.lastIndexOf("/") == -1){
			temp.remove(0,_fileName.lastIndexOf("\\")+1);
			_fileName.remove(_fileName.lastIndexOf("\\")+1, _fileName.size());
		}
		else {
			temp.remove(0,_fileName.lastIndexOf("/")+1);
			_fileName.remove(_fileName.lastIndexOf("/")+1, _fileName.size());
		}
		
		writeFileXml << "	<resdir>" << _fileName.toAscii().data() << "</resdir>" << endl;
		temp.replace(QString(".stl"), QString("_surface.cdb"));
		writeFileXml << "	<res_surface_mesh>" << temp.toAscii().data() << "</res_surface_mesh>" << endl;
		temp.replace(QString("_surface.cdb"), QString("_volume.cdb"));
		writeFileXml << "	<res_volume_mesh>" << temp.toAscii().data() << "</res_volume_mesh>" << endl;
		
		writeFileXml << "</morphing>" << endl;
		outFileXml.close();

	}
}