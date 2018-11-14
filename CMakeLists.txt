cmake_minimum_required(VERSION 2.6)
PROJECT(tools)

FIND_PACKAGE(ITK)
INCLUDE(${USE_ITK_FILE})

FIND_PACKAGE(VTK)
INCLUDE(${VTK_USE_FILE})

# use what QVTK built with
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
  ${tools_SOURCE_DIR}/include
  ${tools_SOURCE_DIR}/include/utils
  ${QT_INCLUDE_DIR}
  ${QT_QTSQL_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  $ENV{INSIGHT_APPLICATIONS}/Auxiliary/vtk
 )
  
SET(ProjectSrcs
  ${tools_SOURCE_DIR}/src/Main.cxx
  ${tools_SOURCE_DIR}/src/MainWindow.cxx
  ${tools_SOURCE_DIR}/src/MeshBoneImplantFitting.cxx
  ${tools_SOURCE_DIR}/src/MeshReader.cxx
  ${tools_SOURCE_DIR}/src/MeshReaderAbaqus.cxx
  ${tools_SOURCE_DIR}/src/MeshReaderMorpherSurface.cxx
  ${tools_SOURCE_DIR}/src/MeshReaderMorpherVolume.cxx
  ${tools_SOURCE_DIR}/src/MeshReaderNeutralFormat.cxx
  ${tools_SOURCE_DIR}/src/MeshWriter.cxx
  ${tools_SOURCE_DIR}/src/MeshWriterAbaqus.cxx
  ${tools_SOURCE_DIR}/src/MeshWriterAnsys.cxx
  ${tools_SOURCE_DIR}/src/MeshSimplifyAndSmooth.cxx
  ${tools_SOURCE_DIR}/src/MeshExtractOuterSurface.cxx  
  ${tools_SOURCE_DIR}/src/MeshMorpherPointPropagation.cxx
  ${tools_SOURCE_DIR}/src/MeshQuality.cxx  
  ${tools_SOURCE_DIR}/src/PointReader.cxx
  ${tools_SOURCE_DIR}/src/PointReaderXyz.cxx
  ${tools_SOURCE_DIR}/src/PointReaderXyzId.cxx
  ${tools_SOURCE_DIR}/src/PointWriter.cxx
  ${tools_SOURCE_DIR}/src/PointWriterXyz.cxx
  ${tools_SOURCE_DIR}/src/PointWriterXyzId.cxx
  ${tools_SOURCE_DIR}/src/PointWriterXyzVector.cxx
  ${tools_SOURCE_DIR}/src/PointWriterMorpherXml.cxx
  ${tools_SOURCE_DIR}/src/ImageHandler.cxx 
  ${tools_SOURCE_DIR}/src/ImageHandlerFloat.cxx  
  ${tools_SOURCE_DIR}/src/VectorImageHandler.cxx
  ${tools_SOURCE_DIR}/src/RegistrationProcrustesAlignment.cxx
  ${tools_SOURCE_DIR}/src/RegistrationValidationVolume.cxx
  ${tools_SOURCE_DIR}/src/RegistrationValidationVolumeImage.cxx
  ${tools_SOURCE_DIR}/src/RegistrationValidationVolumeMesh.cxx
  ${tools_SOURCE_DIR}/src/RegistrationValidationSurface.cxx
  ${tools_SOURCE_DIR}/src/RegistrationValidationSurfaceImage.cxx
  ${tools_SOURCE_DIR}/src/RegistrationValidationSurfaceMesh.cxx
  ${tools_SOURCE_DIR}/src/PCA.cxx
  ${tools_SOURCE_DIR}/src/PCAMesh.cxx
  ${tools_SOURCE_DIR}/src/PCAImages.cxx
  ${tools_SOURCE_DIR}/src/FemAbaqusInpWriter.cxx  
  ${tools_SOURCE_DIR}/src/FemAssigner.cxx
  ${tools_SOURCE_DIR}/src/FemAssignerElements.cxx
  ${tools_SOURCE_DIR}/src/FemAssignerNodes.cxx
  ${tools_SOURCE_DIR}/src/FemForce.cxx
  ${tools_SOURCE_DIR}/src/Rendering.cxx
  ${tools_SOURCE_DIR}/src/RenderingCoordinateSystem.cxx
  ${tools_SOURCE_DIR}/src/RenderingOBB.cxx
  ${tools_SOURCE_DIR}/src/RenderingPoint.cxx
  ${tools_SOURCE_DIR}/src/RenderingPointWithColorbar.cxx
  ${tools_SOURCE_DIR}/src/RenderingPointPicker.cxx
  ${tools_SOURCE_DIR}/src/RenderingMesh.cxx
  ${tools_SOURCE_DIR}/src/RenderingQuadraticTetraMesh.cxx
  ${tools_SOURCE_DIR}/src/VnlReader.cxx
  ${tools_SOURCE_DIR}/src/VnlReaderEValues.cxx
  ${tools_SOURCE_DIR}/src/VnlReaderMatrix.cxx
  ${tools_SOURCE_DIR}/src/VnlReaderVector.cxx 
  ${tools_SOURCE_DIR}/src/VnlWriter.cxx
  ${tools_SOURCE_DIR}/src/VnlWriterEValues.cxx
  ${tools_SOURCE_DIR}/src/VnlWriterMatrix.cxx
  ${tools_SOURCE_DIR}/src/VnlWriterVector.cxx
  ${tools_SOURCE_DIR}/src/VnlWriterVector.cxx
  ${tools_SOURCE_DIR}/src/StatisticsBasics.cxx
  ${tools_SOURCE_DIR}/src/StatisticsDistanceCalculator.cxx
  ${tools_SOURCE_DIR}/src/StatisticsEntropy.cxx
)

SET(UIS 
  ${tools_SOURCE_DIR}/ui/MainWindow.ui)

SET(ProjectHeaders  
  ${tools_SOURCE_DIR}/include/itkLabelOverlapMeasuresImageFilter.h
  ${tools_SOURCE_DIR}/include/itkLabelOverlapMeasuresImageFilter.txx
  ${tools_SOURCE_DIR}/include/itkFieldPCAShapeModelEstimator.h
  ${tools_SOURCE_DIR}/include/itkFieldPCAShapeModelEstimator.txx
  ${tools_SOURCE_DIR}/include/MainWindow.h
  ${tools_SOURCE_DIR}/include/MeshBoneImplantFitting.h
  ${tools_SOURCE_DIR}/include/MeshReader.h
  ${tools_SOURCE_DIR}/include/MeshReaderAbaqus.h
  ${tools_SOURCE_DIR}/include/MeshReaderMorpherSurface.h
  ${tools_SOURCE_DIR}/include/MeshReaderMorpherVolume.h
  ${tools_SOURCE_DIR}/include/MeshReaderNeutralFormat.h  
  ${tools_SOURCE_DIR}/include/MeshWriter.h
  ${tools_SOURCE_DIR}/include/MeshWriterAbaqus.h
  ${tools_SOURCE_DIR}/include/MeshWriterAnsys.h
  ${tools_SOURCE_DIR}/include/MeshSimplifyAndSmooth.h
  ${tools_SOURCE_DIR}/include/MeshExtractOuterSurface.h  
  ${tools_SOURCE_DIR}/include/MeshMorpherPointPropagation.h
  ${tools_SOURCE_DIR}/include/MeshQuality.h
  ${tools_SOURCE_DIR}/include/PointReader.h
  ${tools_SOURCE_DIR}/include/PointReaderXyz.h
  ${tools_SOURCE_DIR}/include/PointReaderXyzId.h
  ${tools_SOURCE_DIR}/include/PointWriter.h
  ${tools_SOURCE_DIR}/include/PointWriterXyz.h
  ${tools_SOURCE_DIR}/include/PointWriterXyzId.h
  ${tools_SOURCE_DIR}/include/PointWriterXyzVector.h
  ${tools_SOURCE_DIR}/include/PointWriterMorpherXml.h
  ${tools_SOURCE_DIR}/include/ImageHandler.h
  ${tools_SOURCE_DIR}/include/ImageHandlerFloat.h
  ${tools_SOURCE_DIR}/include/VectorImageHandler.h
  ${tools_SOURCE_DIR}/include/RegistrationProcrustesAlignment.h
  ${tools_SOURCE_DIR}/include/RegistrationValidationVolume.h
  ${tools_SOURCE_DIR}/include/RegistrationValidationVolumeImage.h
  ${tools_SOURCE_DIR}/include/RegistrationValidationVolumeMesh.h
  ${tools_SOURCE_DIR}/include/RegistrationValidationSurface.h
  ${tools_SOURCE_DIR}/include/RegistrationValidationSurfaceImage.h
  ${tools_SOURCE_DIR}/include/RegistrationValidationSurfaceMesh.h
  ${tools_SOURCE_DIR}/include/PCA.h
  ${tools_SOURCE_DIR}/include/PCAMesh.h
  ${tools_SOURCE_DIR}/include/PCAImages.h
  ${tools_SOURCE_DIR}/include/FemAbaqusInpWriter.h
  ${tools_SOURCE_DIR}/include/FemAssigner.h
  ${tools_SOURCE_DIR}/include/FemAssignerElements.h
  ${tools_SOURCE_DIR}/include/FemAssignerNodes.h
  ${tools_SOURCE_DIR}/include/FemForce.h
  ${tools_SOURCE_DIR}/include/Rendering.h
  ${tools_SOURCE_DIR}/include/RenderingOBB.h
  ${tools_SOURCE_DIR}/include/RenderingCoordinateSystem.h
  ${tools_SOURCE_DIR}/include/RenderingPoint.h
  ${tools_SOURCE_DIR}/include/RenderingPointWithColorbar.h
  ${tools_SOURCE_DIR}/include/RenderingPointPicker.h
  ${tools_SOURCE_DIR}/include/RenderingMesh.h
  ${tools_SOURCE_DIR}/include/RenderingQuadraticTetraMesh.h
  ${tools_SOURCE_DIR}/include/VnlReader.h
  ${tools_SOURCE_DIR}/include/VnlReaderEValues.h
  ${tools_SOURCE_DIR}/include/VnlReaderMatrix.h
  ${tools_SOURCE_DIR}/include/VnlReaderVector.h
  ${tools_SOURCE_DIR}/include/VnlWriter.h
  ${tools_SOURCE_DIR}/include/VnlWriterEValues.h
  ${tools_SOURCE_DIR}/include/VnlWriterMatrix.h
  ${tools_SOURCE_DIR}/include/VnlWriterVector.h
  ${tools_SOURCE_DIR}/include/StatisticsBasics.h
  ${tools_SOURCE_DIR}/include/StatisticsDistanceCalculator.h
  ${tools_SOURCE_DIR}/include/StatisticsEntropy.h
)

QT4_WRAP_UI(UIHeaders ${UIS})
QT4_WRAP_CPP(MOCSrcs ${ProjectHeaders} )
ADD_DEFINITIONS(-DQT_GUI_LIBS -DQT_CORE_LIB -DQT3_SUPPORT)
SET_SOURCE_FILES_PROPERTIES(${ProjectSrcs} PROPERTIES
                            OBJECT_DEPENDS "${UIHeaders}")

# It's nice to have the ui in the windows project file...just double click on it
# and designer comes up on that ui file :)
IF (${CMAKE_BUILD_TOOL} MATCHES "msdev")
  SET (ProjectSrcs ${ProjectSrcs} ${UIS})
ENDIF (${CMAKE_BUILD_TOOL} MATCHES "msdev")
IF (${CMAKE_BUILD_TOOL} MATCHES "devenv")
  SET (ProjectSrcs ${ProjectSrcs} ${UIS})
ENDIF (${CMAKE_BUILD_TOOL} MATCHES "devenv")

ADD_EXECUTABLE( tools MACOSX_BUNDLE ${ProjectHeaders} ${ProjectSrcs} ${UISrcs} ${MOCSrcs})

TARGET_LINK_LIBRARIES( tools
  QVTK
  ${QT_LIBRARIES}
  ${QT_QTSQL_LIBRARIES}
  vtkRendering
  vtkGraphics
  vtkIO
  vtkCommon
  vtkFiltering
  vtkHybrid
  vtkWidgets
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKNumerics
  ITKIOReview
  itkStatistics
)
