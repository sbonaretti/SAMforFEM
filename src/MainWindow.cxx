/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MainWindow.h>

#include <time.h> 

#include <MeshBoneImplantFitting.h>
#include <MeshExtractOuterSurface.h>
#include <MeshQuality.h>
#include <MeshMorpherPointPropagation.h>
#include <MeshReaderAbaqus.h>
#include <MeshReaderMorpherSurface.h>
#include <MeshReaderMorpherVolume.h>
#include <MeshReaderNeutralFormat.h>
#include <MeshSimplifyAndSmooth.h>
#include <MeshWriterAbaqus.h>
#include <MeshWriterAnsys.h>

#include <PointReaderXyz.h>
#include <PointReaderXyzId.h>
#include <PointWriterXyz.h>
#include <PointWriterXyzId.h>
#include <PointWriterMorpherXml.h>

#include <ImageHandler.h>
#include <ImageHandlerFloat.h>
#include <VectorImageHandler.h>

#include <RegistrationProcrustesAlignment.h>
#include <RegistrationValidationSurfaceImage.h>
#include <RegistrationValidationSurfaceMesh.h>
#include <RegistrationValidationVolumeImage.h>
#include <RegistrationValidationVolumeMesh.h>

#include <PCAMesh.h>
#include <PCAImages.h>

#include <FemAbaqusInpWriter.h>
#include <FemAssignerElements.h>
#include <FemAssignerNodes.h>
#include <FemForce.h>

#include <RenderingCoordinateSystem.h>
#include <RenderingOBB.h>
#include <RenderingMesh.h>
#include <RenderingPoint.h>
#include <RenderingPointWithColorbar.h>
#include <RenderingPointPicker.h>
#include <RenderingQuadraticTetraMesh.h>

#include <StatisticsBasics.h>
#include <StatisticsDistanceCalculator.h>

#include <VnlReaderVector.h>
#include <VnlWriterVector.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkLineSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>

#include <QDirModel>
#include <QFile>
#include <QFiledialog>
#include <QModelIndex>
#include <QProcess>
#include <QString>
#include <QStringList>
#include <QTextStream>

#include <vnl/vnl_vector.h>
#include <VnlWriterMatrix.h>//
#include <vnl/vnl_matrix.h>//
#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>


using namespace std;
using namespace image;
using namespace mesh;
using namespace points;
using namespace registration;
using namespace pca;
using namespace fem;
using namespace rendering;
using namespace vnl;
using namespace statistics;


// Constructor
MainWindow::MainWindow(QWidget* parent) 
 : QMainWindow(parent)
{
	// set up gui
	setupUi(this);
	
	
	// data members initialization
	_nOfImageInstances = 0;
	_nOfImageModes = 0;
	_simplifyValue = 0.40;
	_smoothNofIteration = 150;
	_modelVariability = 0.95;
	_nOfMeshInstances = 0;
	_nOfMeshModes = 0;
	_nOfInstances = 0;//
	_modeInterval = 0;//
	_flagMeshEigenvalue = false;
	_flagMeshEigenvector = false;
	_flagMeshAverage = false;
	_PCAflag = 0; // initialized to shape pca
	_assignmentLawOne = 6.850;
	_assignmentLawTwo = 1.49;
	_nodeElementFlag = 0;
	_fileType = 0;
	_positionFlag = 0;
	_loadMagnitude = 0.0;
	_BCtype = 0;
	_meshFlag = false;
	_mechPropFlag = false;
	_loadFlag = false;
	_bCInp = false;

	
	// picker
	_pointPicker = new RenderingPointPicker;
	_cellPicker = vtkCellPicker::New();
    _cellPicker->SetTolerance(0.01);
	qvtkWidget->GetInteractor()->SetPicker(_cellPicker);

	
	// Image tool - slots
	connect(ImageResample, SIGNAL(clicked()), this, SLOT(resample()));
	connect(ItkToVtkToItk, SIGNAL(clicked()), this, SLOT(itkToVtkToItk()));
	connect(PseudoCalibration, SIGNAL(clicked()), this, SLOT(pseudoCalibration()));
	connect(MaskGreyLevel, SIGNAL(clicked()), this, SLOT(maskGreyLevel()));
	connect(MaskImage, SIGNAL(clicked()), this, SLOT(maskImage()));
	connect(BinaryMaskDcm, SIGNAL(clicked()), this, SLOT(binaryMaskDcm()));

	connect(FindReference, SIGNAL(clicked()), this, SLOT(findReference()));
	connect(VfToDvf, SIGNAL(clicked()), this, SLOT(vfToDvf()));
	connect(RecreateMoving, SIGNAL(clicked()), this, SLOT(recreateMoving()));
	connect(VfToDvfInverted, SIGNAL(clicked()), this, SLOT(vfToDvfInverted()));
	connect(RecreateReference, SIGNAL(clicked()), this, SLOT(recreateReference()));

	connect(SurfaceComparisonImage, SIGNAL(clicked()), this, SLOT(surfaceComparisonImage()));
	connect(VolumeComparisonImage, SIGNAL(clicked()), this, SLOT(volumeComparisonImage()));
	
	connect(ShapeImagePCA, SIGNAL(toggled(bool)), this, SLOT(shapeImagePCA()));
	connect(IntensityImagePCA, SIGNAL(toggled(bool)), this, SLOT(intensityImagePCA()));
	connect(CombinedImagePCA, SIGNAL(toggled(bool)), this, SLOT(combinedImagePCA()));
	connect(CalculateImagePCA, SIGNAL(clicked()), this, SLOT(calculateImagePCA()));

	connect(RecreateImageInstance, SIGNAL(clicked()), this, SLOT(recreateImageInstance()));
		
	connect(NOfImageInstances, SIGNAL(valueChanged(int)), this, SLOT(nOfImageInstances(int)));
	connect(NOfImageModes_2, SIGNAL(valueChanged(int)), this, SLOT(nOfImageModes(int)));
	connect(CreateImageInstances, SIGNAL(clicked()), this, SLOT(createImageInstances()));

	
	// Mesh tool - slots
	connect(LoadImageMesh, SIGNAL(clicked()), this, SLOT(loadImageMesh()));
	connect(Simplify, SIGNAL(valueChanged(double)), this, SLOT(simplify(double)));
	connect(Smooth, SIGNAL(valueChanged(int)), this, SLOT(smooth(int)));
	connect(CreateMesh, SIGNAL(clicked()), this, SLOT(createMesh()));
	
	connect(LoadSurfaceMesh, SIGNAL(clicked()), this, SLOT(loadSurfaceMesh())); 
	connect(CreateVolumeMesh, SIGNAL(clicked()), this, SLOT(createVolumeMesh()));
	
	connect(MeshJacobian, SIGNAL(clicked()), this, SLOT(meshJacobian()));
	connect(MeshEdgeRatio, SIGNAL(clicked()), this, SLOT(meshEdgeRatio()));
	connect(MeshAngle, SIGNAL(clicked()), this, SLOT(meshAngle()));
	
	connect(PickLandmarks, SIGNAL(clicked()), this, SLOT(pickLandmarks()));
	connect(CreateXml, SIGNAL(clicked()), this, SLOT(createXml()));
	connect(CreateXmlFromFiles, SIGNAL(clicked()), this, SLOT(createXmlFromFiles()));
	connect(RunMorpher, SIGNAL(clicked()), this, SLOT(runMorpher()));
	connect(SurfaceValidation, SIGNAL(clicked()), this, SLOT(surfaceValidation()));
	connect(VolumeValidation, SIGNAL(clicked()), this, SLOT(volumeValidation()));
	
	connect(ExtractCoordinatesAndIntensities, SIGNAL(clicked()), this, SLOT(extractCoordinatesAndIntensities()));
	connect(ExtractNodesAndElements, SIGNAL(clicked()), this, SLOT(extractNodesAndElements()));
	connect(MeshAlignment, SIGNAL(clicked()), this, SLOT(meshAlignment()));
	

	connect(ShapeMeshPCA, SIGNAL(toggled(bool)), this, SLOT(shapeMeshPCA()));
	connect(IntensityMeshPCA, SIGNAL(toggled(bool)), this, SLOT(instensityMeshPCA()));
	connect(CombinedMeshPCA, SIGNAL(toggled(bool)), this, SLOT(combinedMeshPCA()));
	connect(CalculateMeshPCA, SIGNAL(clicked()), this, SLOT(calculateMeshPCA()));
	connect(RecreateMeshInstance, SIGNAL(clicked()), this, SLOT(recreateMeshInstance()));
	connect(NOfMeshInstances, SIGNAL(valueChanged(int)), this, SLOT(nOfMeshInstances(int)));
	connect(NOfMeshModes, SIGNAL(valueChanged(int)), this, SLOT(nOfMeshModes(int))); 
	connect(CreateMeshInstances, SIGNAL(clicked()), this, SLOT(createMeshInstances()));

	// FEM tool - slots
	connect(LoadMesh, SIGNAL(clicked()), this, SLOT(loadMesh()));
	connect(LoadImageAssigner, SIGNAL(clicked()), this, SLOT(loadImage()));
	connect(AssignmentLawFirst, SIGNAL(valueChanged(double)), this, SLOT(assignmentLawFirst(double)));
	connect(AssignmentLawSecond, SIGNAL(valueChanged(double)), this, SLOT(assignmentLawSecond(double)));
	connect(AssignToElements, SIGNAL(toggled(bool)), this, SLOT(assignToElements()));
	connect(AssignToNodes, SIGNAL(toggled(bool)), this, SLOT(assignToNodes()));
	connect(AbaqusMatProp, SIGNAL(toggled(bool)), this, SLOT(abaqusMatProp()));
	connect(TxtMatProp, SIGNAL(toggled(bool)), this, SLOT(txtMatProp()));
	connect(Assign, SIGNAL(clicked()), this, SLOT(assign()));

	connect(LoadBoneMesh, SIGNAL(clicked()), this, SLOT(loadBoneMesh()));
	connect(LoadImplantMesh, SIGNAL(clicked()), this, SLOT(loadImplantMesh()));
	connect(Proximal, SIGNAL(toggled(bool)), this, SLOT(proximal()));
	connect(Diaphyseal, SIGNAL(toggled(bool)), this, SLOT(diaphyseal()));
	connect(Distal, SIGNAL(toggled(bool)), this, SLOT(distal()));
	connect(ComputeAndSave, SIGNAL(clicked()), this, SLOT(computeAndSave()));

	connect(PickApplicationPoints, SIGNAL(clicked()), this, SLOT(pickApplicationPoints()));
	connect(SaveApplicationPoints, SIGNAL(clicked()), this, SLOT(saveApplicationPoints()));
	connect(PropagateApplicationPoints, SIGNAL(clicked()), this, SLOT(propagateApplicationPoints()));
	connect(Walking, SIGNAL(clicked()), this, SLOT(walking()));
	connect(Falling, SIGNAL(clicked()), this, SLOT(falling()));
	
	
	connect(CreateAbaqusInp, SIGNAL(clicked()), this, SLOT(createAbaqusInp()));
	connect(CreateAbaqusInpFalling, SIGNAL(clicked()), this, SLOT(createAbaqusInpFalling()));
	connect(RunSimulation, SIGNAL(clicked()), this, SLOT(runSimulation()));

	connect(ExtractNeckIntensities, SIGNAL(clicked()), this, SLOT(extractNeckIntensities()));

	
	// Rendering
	connect(RenderMetafile, SIGNAL(triggered()), this, SLOT(renderMetafile()));
	connect(RenderStl, SIGNAL(triggered()), this, SLOT(renderStl()));
	connect(RenderSurfaceCdb, SIGNAL(triggered()), this, SLOT(renderSurfaceCdb()));
	connect(RenderVolumeCdb, SIGNAL(triggered()), this, SLOT(renderVolumeCdb()));
	connect(RenderVolumeInpAbaqus, SIGNAL(triggered()), this, SLOT(renderVolumeInpAbaqus()));
	connect(RenderVtkPolyData, SIGNAL(triggered()), this, SLOT(renderVtkPolyData()));
	connect(RenderPointsXyz, SIGNAL(triggered()), this, SLOT(renderPointsXyz()));
	connect(RenderPointsXyzWithColorbar, SIGNAL(triggered()), this, SLOT(renderPointsXyzWithColorbar()));
	connect(RenderPointsXyzWithColorbarC, SIGNAL(triggered()), this, SLOT(renderPointsXyzWithColorbarColors()));
	connect(RenderLine, SIGNAL(triggered()), this, SLOT(renderLine()));
	connect(RenderOBB, SIGNAL(triggered()), this, SLOT(renderOBB()));
	connect(RenderCoordinateSystem, SIGNAL(triggered()), this, SLOT(renderCoordinateSystem()));
	connect(EuclideanDistancePoints, SIGNAL(triggered()), this, SLOT(euclideanDistancePoints()));
	connect(DrawLineStart, SIGNAL(triggered()), this, SLOT(pickApplicationPoints()));
	connect(DrawLineStop, SIGNAL(triggered()), this, SLOT(drawLineStop()));
	connect(ResetView, SIGNAL(triggered()), this, SLOT(resetView()));


	// Format conversion
	connect(DcmToMhd, SIGNAL(triggered()), this, SLOT(dicomToMetafile()));
	connect(MhdToDcm, SIGNAL(triggered()), this, SLOT(metafileToDicom()));
	connect(MhdToHdr, SIGNAL(triggered()), this, SLOT(metafileToAnalyse()));
	
	connect(MorpherSurfaceCdbToStl_2, SIGNAL(triggered()), this, SLOT(morpherSurfaceCdbToStl()));
	connect(MorpherVolumeCdbToVtk_2, SIGNAL(triggered()), this, SLOT(morpherVolumeCdbToVtk()));
	connect(AbaqusToCdb, SIGNAL(triggered()), this, SLOT(abaqusToCdb()));
	connect(CdbToInpAbaqus, SIGNAL(triggered()), this, SLOT(cdbToInpAbaqus()));
	
	connect(PointsStlToTxt, SIGNAL(triggered()), this, SLOT(ptStlToTxt()));
	connect(PointsVtkToTxt, SIGNAL(triggered()), this, SLOT(ptVtkToTxt()));
	
	
	// QT/VTK interact
	_ren = vtkRenderer::New();
	_background[0] = 1.0; _background[1] = 1.0; _background[2] = 1.0;
	_ren->SetBackground(_background); // Background color white
	qvtkWidget->GetRenderWindow()->AddRenderer(_ren);
	
	
	// Reset camera
	_ren->ResetCamera();
	qvtkWidget->GetRenderWindow()->Render();
	 
};

// Destructor
MainWindow::~MainWindow() {
};

/************************* IMAGE TOOL *************************/
// ad hoc preprocessing
void MainWindow::resample(){

	/*
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load images", "C:/0.Data/test data/image extrusion", "Images (*.mhd)");
		
	for (int i=0; i<fileNames.size(); i++){
	
		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImageFileName(fileNames[i]);
		imageHandler->MetafileReaderUpdate();
		imageHandler->ResampleUpdate();
		fileNames[i].replace(QString(".mhd"), QString("_resampled.mhd"));

		imageHandler->SetImageFileName(fileNames[i]);
		std::cout << fileNames[i].toAscii().data() << std::endl;
		imageHandler->MetafileWriterUpdate();
		
		std::cout << "resampled image written" << std::endl;

		delete imageHandler;
	}
	*/
	
	/*ImageHandler* test = new ImageHandler;
	test->ExtractMasks();
	*/
	/*ImageHandler* test = new ImageHandler;
	test->CreateTestImages();
	*/
	/*
	// mask combination
	QStringList maskFileNames = QFileDialog::getOpenFileNames(this, "Load images", "C:/0.Data/0. Original data/3. femur left masks/1. mine", "Images (*.mhd)");
	QStringList maskFileNamesTwo = QFileDialog::getOpenFileNames(this, "Load images", "C:/0.Data/0. Original data/3. femur left masks/2. stryker modified", "Images (*.mhd)");
	
	QString folderName = QFileDialog::getExistingDirectory(this, "Select parent folder", "C:/0.Data/0. Original data/0. femur left dicom/0. femur left"); 
	
	QStringList list;
	QDir dir(folderName);
	QStringList filters;
	filters << "*Segmented-SLICES";
	dir.setNameFilters(filters);
	list = dir.entryList();

	QStringList segmentedNames;
	segmentedNames = list; // just to get the size
	std::cout << segmentedNames.size() << std::endl;

	for (int i=0; i<list.size(); i++){
		QString temp = folderName;
		temp.append("\\");
		segmentedNames[i] = temp.append(list[i]); 
	}

	ImageHandler* imageHandler = new ImageHandler;
	imageHandler->SetImageFileNames(segmentedNames);
	imageHandler->SetMaskFileNames(maskFileNames);
	imageHandler->SetMaskFileNamesTwo(maskFileNamesTwo);
	imageHandler->CombineMasks();
	*/
	
	
	// histogram 
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load images", "C:/0.Data/1. Registered data/1a. dvf and warped to reference/images", "Images (*.mhd)");
		 
	ImageHandler* imageHandler = new ImageHandler;
	imageHandler->SetImageFileNames(fileNames);
	imageHandler->ImageHistogram();
		
	delete imageHandler;
	
	
	/*
	// cropping
	QString fileName = QFileDialog::getOpenFileName(this, "Load image", " ", "Image (*.mhd)");
	ImageHandler* imageHandler = new ImageHandler;
	imageHandler->SetImageFileName(fileName);
	imageHandler->MetafileReaderUpdate();
	imageHandler->Crop();
	std::cout << "done" << std::endl;
	*/

	/*
	QString fileName = ("C:/0.Data/test data/boxplot/randomData2.txt");
	vnl_vector<double> test;
	
	VnlReaderVector* reader = new VnlReaderVector;
	reader->SetFileName(fileName);
	reader->Update();
	test = reader->GetVnlVector();
	
	StatisticsBasics* calc = new StatisticsBasics;
	calc->SetVector(test);
	calc->CalculateBoxPlot();
	*/

	/*
	// dicom header reader
	ImageHandler* imageHandler = new ImageHandler; 
	QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"), 
		"/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
	if(dir.compare("") == 0) return;
		
	// folders in the dir
	QDirModel *A = new QDirModel;
    QModelIndex indexA = A->index(dir);
	int numB = A->rowCount(indexA);
	
	QDirModel *B = new QDirModel;
	
    for (int i = 0; i < numB; ++i) {
			
		// subfolder in the folder 
		QModelIndex indexB = A->index(i, 0, indexA);
		QString folder = A->data(indexB, Qt::DisplayRole).toString(); // inner folder name
		std::cout << "\n" << i+1 << " " << folder.ascii() << std::endl;		
		QString totalpath =  dir + "/" + folder; // complete path of each folder
		QModelIndex indexC = B->index(totalpath);
		int numFiles = B->rowCount(indexC);
		
		for (int j=0; j<numFiles; j++){
			QModelIndex fileIndex = B->index(j, 0, indexC);
			QString fileName = B->data(fileIndex, Qt::DisplayRole).toString(); // file name
			if (fileName.endsWith(".dcm")) {
				QString totalpathFile = totalpath + "/" + fileName; // complete file name
				imageHandler->SetImageFileName(totalpathFile.ascii());
				imageHandler->DicomHeaderReaderFromSlice();
				break;
			}
			else j++;
		}
	}
	*/

	/*
	// subtract images
	// load first image
	QString fileName = QFileDialog::getOpenFileName(this, "Load image", " ", "Image (*.mhd)");
	ImageHandler* imageHandler = new ImageHandler;
	imageHandler->SetImageFileName(fileName);
	std::cout << fileName.toAscii().data() << std::endl;
	imageHandler->MetafileReaderUpdate();

	// load second image
	QString fileName2 = QFileDialog::getOpenFileName(this, "Load image", " ", "Image (*.mhd)");
	ImageHandler* imageHandler2 = new ImageHandler;
	imageHandler2->SetImageFileName(fileName2);
	std::cout << fileName2.toAscii().data() << std::endl;
	imageHandler2->MetafileReaderUpdate();
	
	// subtract images
	ImageHandler* imageHandler3 = new ImageHandler;
	imageHandler3->SetImage(imageHandler->GetImage());
	imageHandler3->SetImage2(imageHandler2->GetImage());
	imageHandler3->SubtractImages();
	std::cout << "subtracted image written" << std::endl;
	*/

	
}
void MainWindow::itkToVtkToItk(){

	// to improve
	QString fileName = QFileDialog::getOpenFileName(this, "Load image", "*.mhd");
	
	ImageHandler* imageHandler = new ImageHandler;
	imageHandler->SetImageFileName(fileName);
	imageHandler->MetafileReaderUpdate();
	fileName.replace(QString(".mhd"), QString("_vtkCoor.mhd"));
	imageHandler->SetImageFileName(fileName);
	imageHandler->ITKtoVTKtoITK();

}
void MainWindow::pseudoCalibration(){

	std::cout << "! all images will be multiplied by 1000 in order to have short images - from [mg/mm3 to g/mm3]" << std::endl;

	/*
	// min selected from the marrow; max selected from the image
	// the reference masked is used to find the minimum of the marrow area of each bone

	// load images
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load images", "C:/0.Data/test data/fracture risk/data", "Images (*.mhd)");
	// load the correspondent svf 
	QStringList svfFileNames = QFileDialog::getOpenFileNames(this, "Load svf","C:/0.Data/test data/fracture risk/data/0.referenceDVF","Images(*.mhd)");
	// load reference mask
	QString referenceMaskFileName = QFileDialog::getOpenFileName(this, "Load reference mask","C:/0.Data/test data/fracture risk/data/0.referenceDVF/mask","Images(*.mhd)");
	
	// load the reference mask
	ImageHandlerFloat* maskHandler = new ImageHandlerFloat;
	maskHandler->SetImageFileName(referenceMaskFileName);
	maskHandler->MetafileReaderUpdate();


	for (int i=0; i<fileNames.size(); i++){
		
		// load image to calibrate
		ImageHandlerFloat* imageHandlerFloat = new ImageHandlerFloat;
		std::cout << fileNames[i].toAscii().data() << std::endl;
		imageHandlerFloat->SetImageFileName(fileNames[i]);
		imageHandlerFloat->MetafileReaderUpdate();

		// load the correspondent svf
		VectorImageHandler* vectorImageHandler = new VectorImageHandler;
		QString temp = fileNames[i];
		temp.remove(0, temp.lastIndexOf("_"));
		temp.remove(temp.size()-3, temp.size()-1);

		for (int y=0; y<svfFileNames.size(); y++){
			if (svfFileNames[y].contains(temp)){
				vectorImageHandler->SetFieldFileName(svfFileNames[y]);
				std::cout << svfFileNames[y].ascii() << std::endl;
			}
		}
		vectorImageHandler->MetafileReaderUpdate();
		
		std::cout << "from SVT to DVF" << std::endl;
		vectorImageHandler->VFtoDVFnotInverted();
		
		imageHandlerFloat->SetSVF(vectorImageHandler->GetField());

		// set the reference mask
		imageHandlerFloat->SetMask(maskHandler->GetImage());
		
		// find min
		std::cout << "find marrow minimum" << std::endl;
		//imageHandlerFloat->FindMinMarrow();
		
		// calibrate
		imageHandlerFloat->PseudoCalibration();

		// write image
		fileNames[i].replace(QString(".mhd"), QString("_pseudoCalibrated.mhd"));
		imageHandlerFloat->SetImageFileName(fileNames[i]);
		std::cout << fileNames[i].toAscii().data() << std::endl;
		imageHandlerFloat->MetafileWriterUpdate();
		
		std::cout << "pseudo-calibrated image written" << std::endl;
		std::cout << std::endl;

		delete imageHandlerFloat;
		delete vectorImageHandler;

	}

	delete maskHandler;
	*/



	
	// load images
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load images", "C:/0.Data/2011_FemaleMale_Study/6. new instances/male/2.images", "Images (*.mhd)");
	
	// if no reference mask is used, set the min in other way in ImageHandlerFloat::PseudoCalibration
	for (int i=0; i<fileNames.size(); i++){
	
		ImageHandlerFloat* imageHandler = new ImageHandlerFloat;
		imageHandler->SetImageFileName(fileNames[i]);
		imageHandler->MetafileReaderUpdate();
		imageHandler->PseudoCalibration();

		/* // for originals
		QString temp = fileNames[i];
		if (temp.lastIndexOf("/") == -1){
				//temp.remove("inp");
				temp.remove(0, temp.lastIndexOf("\\")+1);
				temp.replace(QString("."), QString("_"));
		}
		else {
			//temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.replace(QString("."), QString("_"));
		}
		int length = temp.length();
		temp.remove(5,length);
		temp.append("_calibrated.mhd");
		
		imageHandler->SetImageFileName(temp);
		std::cout << temp.toAscii().data() << std::endl;
		imageHandler->MetafileWriterUpdate();
		*/

		// for instances
		fileNames[i].replace(QString(".mhd"), QString("_pseudoCalibrated.mhd"));
		imageHandler->SetImageFileName(fileNames[i]);
		std::cout << fileNames[i].toAscii().data() << std::endl;
		imageHandler->MetafileWriterUpdate();
		
		
		std::cout << "pseudo-calibrated image written" << std::endl;

		delete imageHandler;
	}
	
}
void MainWindow::maskGreyLevel(){

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load the masked images", "C:/0.Data/1. Registered data/1. registered mhd masks - ref 00705/3. 1+2 all registered mhd masks",  "*warpedToRef.mhd");
	 
	for (int i=0; i<fileNames.size(); i++){
	
		std::cout << "mask: " << fileNames[i].ascii() << std::endl;
	
		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImageFileName(fileNames[i]);
		imageHandler->MetafileReaderUpdate();
		imageHandler->MaskGreyLevelsUpdate();
		fileNames[i].replace(QString(".mhd"), QString("_noTransitions.mhd"));
		imageHandler->SetImageFileName(fileNames[i]);
		imageHandler->MetafileWriterUpdate();

		delete imageHandler;


	}

}
void MainWindow::maskImage(){

	QString ctFileName = QFileDialog::getOpenFileName(this, "Load ct image", "*.mhd");
	QString maskFileName = QFileDialog::getOpenFileName(this, "Load mask image", "*.mhd");

	ImageHandler* imageHandler = new ImageHandler;
	imageHandler->SetImageFileName(ctFileName);
	imageHandler->MetafileReaderUpdate();
	imageHandler->SetMaskFileName(maskFileName);
	imageHandler->MaskReaderUpdate();
	imageHandler->MaskToImage();
	ctFileName.replace(QString(".mhd"), QString("_segmented.mhd"));
	imageHandler->SetImageFileName(ctFileName);
	imageHandler->MetafileWriterUpdate();
}
void MainWindow::binaryMaskDcm(){

	
	std::cout << "mask image" << std::endl;
	std::cout << std::endl;

	QString folderName = QFileDialog::getExistingDirectory(this, "Select image folders"); 
	QStringList list;
	QDir dir(folderName);
	QStringList filters;
	filters << "*SLICES";
	dir.setNameFilters(filters);
	list = dir.entryList();
		
	for(int i = 0; i < list.size(); ++i) {

		ImageHandler* imageHandler = new ImageHandler;
		QString folderName;
		
		// image
		folderName = (dir.absolutePath() + "/" + list[i]);
		std::cout << "image: " << folderName.toAscii().data() << std::endl;
		imageHandler->SetDicomSerieFolderName(folderName);
		imageHandler->DicomReaderUpdate();

		// mask
		folderName.replace(QString ("Segmented"), QString("Classified"));
		folderName.append("-ReClassified");
		std::cout << "mask: " << folderName.toAscii().data()<< std::endl;
		imageHandler->SetDicomMaskFolderName(folderName);
		imageHandler->DicomMaskReaderUpdate();

		// update
		imageHandler->MaskToImage();

		// overwrite image
		folderName = (dir.absolutePath() + "/" + "ReSegmented" + "/" + list[i]);
		std::cout << folderName.toAscii().data()<< std::endl;
		dir.mkdir(folderName);
		imageHandler->SetDicomSerieFolderName(folderName);
		imageHandler->DicomWriterUpdate();
		
		delete imageHandler;
	}

	std::cout << "images masked" << std::endl;


	


}

// correspondeces - registration
void MainWindow::findReference(){
	
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load velocity fields", " ", "Velocity fields (*_velocity_field.mhd)");
	
	VectorImageHandler* vectorImageHandler = new VectorImageHandler;
	vectorImageHandler->SetFieldFileNames(fileNames);
	vectorImageHandler->SelectReferenceAsMinimumDistanceToAverage();

}
void MainWindow::vfToDvf(){

	std::cout << "velocity field to DVF" << std::endl;

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load velocity fields", " ", "*.mhd");//"Velocity fields (*_velocity_field.mhd)");
	 
	VectorImageHandler* vectorImageHandler = new VectorImageHandler;
	
	for (int i=0; i<fileNames.size(); i++){
		
		std::cout << fileNames[i].ascii() << std::endl;
		vectorImageHandler->SetFieldFileName(fileNames[i]);
		vectorImageHandler->MetafileReaderUpdate();
		vectorImageHandler->VFtoDVFnotInverted();
		fileNames[i].replace(QString(".mhd"), QString("_dvf.mhd"));
		vectorImageHandler->SetFieldFileName(fileNames[i]);
		vectorImageHandler->MetafileWriterUpdate();
	}

	delete vectorImageHandler;

}
void MainWindow::recreateReference(){

	
	QStringList imageFileNames = QFileDialog::getOpenFileNames(this, "Load images", " ", "*_rigid.mhd");
	
	for (int i=0; i<imageFileNames.size(); i++){

		// read the image 
		std::cout << "image: " << imageFileNames[i].ascii() << std::endl;
		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImageFileName(imageFileNames[i]);
		imageHandler->MetafileReaderUpdate();
		
		// read the DVF
		imageFileNames[i].replace(QString("_rigid.mhd"), QString("_velocity_field_dvf.mhd"));
		VectorImageHandler* vectorHandler = new VectorImageHandler;
		std::cout << "DVF: " << imageFileNames[i].ascii() << std::endl;
		vectorHandler->SetFieldFileName(imageFileNames[i]);
		vectorHandler->MetafileReaderUpdate();
			
		//apply DVF to image
		std::cout << "apply DVF to image" << std::endl;
		imageHandler->SetDeformationVectorField(vectorHandler->GetField());
		imageHandler->ApplyDVFToImageUpdate();

		// save image
		imageFileNames[i].replace(QString("_velocity_field_dvf.mhd"), QString("_rigid_warpedToRef.mhd"));
		std::cout << "output image: " << imageFileNames[i].ascii() << std::endl;
		imageHandler->SetImageFileName(imageFileNames[i]);
		std::cout << imageFileNames[i].ascii() << std::endl;
		imageHandler->MetafileWriterUpdate();

		//cleaning up
		delete vectorHandler;
		delete imageHandler;
	}
	
		
	/**** MASKS WARPED TO REFERENCE ****/
	/*
	// for masks
	QStringList imageFileNames = QFileDialog::getOpenFileNames(this, "Load masks", "C:/0.Data/1. Registered data/1. registered mhd masks - ref 00705/3. 1+2 all registered mhd masks", "*ridig_mask.mhd"); 
	QStringList dvfs = QFileDialog::getOpenFileNames(this, "Load dvfs", "C:/0.Data/1. Registered data/1a. dvf and warped to reference/dvfs", "*.mhd"); 
	
	for (int i=0; i<imageFileNames.size(); i++){

		// read the image 
		std::cout << "image: " << imageFileNames[i].ascii() << std::endl;
		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImageFileName(imageFileNames[i]);
		imageHandler->MetafileReaderUpdate();
		
		VectorImageHandler* vectorHandler = new VectorImageHandler;
		
		QString temp = imageFileNames[i];
		if (temp.lastIndexOf("/") == -1){
			temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
			temp.remove(0, temp.lastIndexOf("\\")+1);
		}
		else {
			temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
			temp.remove(0, temp.lastIndexOf("/")+1);
		}
		for (int y=0; y<dvfs.size(); y++){
			if (dvfs[y].contains(temp)){
				std::cout << "DVF: " << dvfs[y].ascii() << std::endl;
				vectorHandler->SetFieldFileName(dvfs[y]);
		vectorHandler->MetafileReaderUpdate();
			}
		}
				
		//apply DVF to image
		std::cout << "apply DVF to image" << std::endl;
		imageHandler->SetDeformationVectorField(vectorHandler->GetField());
		imageHandler->ApplyDVFToImageUpdate();

		// save image
		imageFileNames[i].replace(QString(".mhd"), QString("-warpedToRef.mhd"));
		std::cout << "output image: " << imageFileNames[i].ascii() << std::endl;
		imageHandler->SetImageFileName(imageFileNames[i]);
		std::cout << imageFileNames[i].ascii() << std::endl;
		imageHandler->MetafileWriterUpdate();

		//cleaning up
		delete vectorHandler;
		delete imageHandler;
	}
	*/
		
}
void MainWindow::vfToDvfInverted(){

	std::cout << "velocity field to inverted DVF" << std::endl;

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load velocity fields", " ", "Velocity fields (*.mhd)");
	 
	VectorImageHandler* vectorImageHandler = new VectorImageHandler;
	
	for (int i=0; i<fileNames.size(); i++){
		
		std::cout << fileNames[i].ascii() << std::endl;
		vectorImageHandler->SetFieldFileName(fileNames[i]);
		vectorImageHandler->MetafileReaderUpdate();
		vectorImageHandler->VFtoDVFinverted();
		fileNames[i].replace(QString(".mhd"), QString("_dvf_inverted.mhd"));
		vectorImageHandler->SetFieldFileName(fileNames[i]);
		vectorImageHandler->MetafileWriterUpdate();
	}

	delete vectorImageHandler;
}
void MainWindow::recreateMoving(){


	QString refFileName = QFileDialog::getOpenFileName(this, "Load reference", " ", "*.mhd");
	QStringList dvfFileNames = QFileDialog::getOpenFileNames(this, "Load deformation vector fields", " ", "*.mhd");
			
	for (int i=0; i<dvfFileNames.size(); i++){

		// read the reference (inside the loop to avoid that the warped image is kept in memory)
		std::cout << "reference: " << refFileName.toAscii().data() << std::endl;
		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImageFileName(refFileName);
		imageHandler->MetafileReaderUpdate();
	
		// read the DVF
		std::cout << "DVF: " << dvfFileNames[i].ascii() << std::endl;
		VectorImageHandler* vectorHandler = new VectorImageHandler;
		vectorHandler->SetFieldFileName(dvfFileNames[i]);
		vectorHandler->MetafileReaderUpdate();
		
		//apply DVF to image
		std::cout << "apply DVF to image" << std::endl;
		imageHandler->SetDeformationVectorField(vectorHandler->GetField());
		imageHandler->ApplyDVFToImageUpdate();

		// save image
		dvfFileNames[i].replace(QString(".mhd"), QString("_rigid_warpedFromRef.mhd"));
		imageHandler->SetImageFileName(dvfFileNames[i]);
		std::cout << dvfFileNames[i].ascii() << std::endl;
		imageHandler->MetafileWriterUpdate();

		//cleaning up
		delete vectorHandler;
		delete imageHandler;
	}

	
	
}
void MainWindow::surfaceComparisonImage(){

	QStringList registeredFileNames = QFileDialog::getOpenFileNames(this, "Load registered", "C:/0.Data/1. Registered data/1a. dvf-1 and warped to originals/marching cubes .stl", "*.stl");
	QStringList originalFileNames = QFileDialog::getOpenFileNames(this, "Load original images", "C:/0.Data/1. Registered data/1. registered mhd - ref 00705/marching cubes .stl", "*.stl");
	
	RegistrationValidationSurfaceMesh* validation = new RegistrationValidationSurfaceMesh;
	validation->SetFlag(1);
	validation->SetRegisteredFileNames(registeredFileNames);
	validation->SetOriginalFileNames(originalFileNames);
	validation->Update();
		
	delete validation;

}
void MainWindow::volumeComparisonImage(){
	
	QString referenceMaskFileName = QFileDialog::getOpenFileName(this, "Load reference mask", "C:/0.Data/1. Registered data/0. ref 00705", "*.mhd");
	QStringList warpedToReferenceMaskFileNames = QFileDialog::getOpenFileNames(this, "Load images warped to reference", "C:/0.Data/1. Registered data/1. registered mhd masks - ref 00705/3. 1+2 all registered mhd masks", "*warpedToRef.mhd");
	
	RegistrationValidationVolumeImage* validation = new RegistrationValidationVolumeImage;
	validation->SetReferenceMaskFileName(referenceMaskFileName);
	validation->SetMovingMaskFileNames(warpedToReferenceMaskFileNames);
	validation->Update();

	delete validation;

}

// statistical model 
void MainWindow::shapeImagePCA(){
	
	_PCAflag = 1;

}
void MainWindow::intensityImagePCA(){

	_PCAflag = 2;

}
void MainWindow::combinedImagePCA(){

	_PCAflag = 3;

}
void MainWindow::calculateImagePCA(){

	
	// shape PCA
	if (_PCAflag == 1){
		
		std::cout << "compute shape PCA" << std::endl;
		
		QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load velocity fields", "C:/0.Data/1. Registered data/1. registered mhd - ref 00705", "Velocity fields (*_velocity_field.mhd);; (*.mhd)");
		
		PCAImages* PCA = new PCAImages;
		PCA->SetImageFileNames(fileNames);
		PCA->ShapePCA();
		
		//delete PCA;
	}
	
	// instensity PCA
	else if (_PCAflag == 2){
				
		std::cout << "! images must be warped to the reference shape" << std::endl;

		QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load images", "C:/0.Data/1. Registered data/1a. dvf and warped to reference/images", "Images (*.mhd)");
		
		PCAImages* PCA = new PCAImages;
		PCA->SetImageFileNames(fileNames);
		PCA->IntensityPCA();
		
		delete PCA;
		
	}

	// combined PCA
	else if (_PCAflag == 3){
		
		std::cout << "compute combined PCA" << std::endl;

		std::cout << "!!! vfs and warped images have to be loaded in the same order as they were for the shape and intensity pca !!!" << std::endl;

		QStringList shapeDatasetFileNames = QFileDialog::getOpenFileNames(this, "Load vf files (same order as for the shape pca)","C:/0.Data/3. Female Male/1. dvf/male", "*.mhd");
		QString shapeAverageFileName = QFileDialog::getOpenFileName(this, "Load shape average file","C:/0.Data/3. Female Male/3. shape pca/male", "*.mhd");
		QStringList shapeEvectorsFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors files","C:/0.Data/3. Female Male/3. shape pca/male", "*.mhd");
		
		QStringList intensityDatasetFileName = QFileDialog::getOpenFileNames(this, "Load warpedToRef files (same order as for the intensity pca)","C:/0.Data/3. Female Male/2. warped to reference/male", "*.mhd");
		QString intensityAverageFileName = QFileDialog::getOpenFileName(this, "Load intensity average file","C:/0.Data/3. Female Male/4. intensity pca/male", "*.mhd");
		QStringList intensityEvectorsFileNames = QFileDialog::getOpenFileNames(this, "Load intensity eigenvectors files","C:/0.Data/3. Female Male/4. intensity pca/male", "*.mhd");
		
		QString wFileName = QFileDialog::getOpenFileName(this, "Load w file name", "C:/0.Data/3. Female Male/5. combined pca/male", "(*.txt)");
		

		double start = clock();

		PCAImages* PCA = new PCAImages;

		PCA->SetVelocityFieldFileNames(shapeDatasetFileNames);
		PCA->SetShapeAverageFileName(shapeAverageFileName);
		PCA->SetShapeEVectorsFileNames(shapeEvectorsFileNames);

		PCA->SetImageFileNames(intensityDatasetFileName);
		PCA->SetIntensityAverageFileName(intensityAverageFileName);
		PCA->SetIntensityEVectorsFileNames(intensityEvectorsFileNames);

		PCA->SetWFileName(wFileName);

		PCA->CombinedPCA();

		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "Computation time: " << total << " sec. (" << total/60 << " min.)" << std::endl;
	
	
	}

	else if (_PCAflag != 1 && _PCAflag != 2 && _PCAflag != 3){

		/*
		typedef float VoxelType;
		static const unsigned int Dimension = 3;
		typedef itk::Image< VoxelType, Dimension > ImageType2;
		typedef itk::ImageFileReader< ImageType2 > ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName("C:/0.Data/test data/debug image-mesh pca/images/bone 1/dvf test.mhd");
		reader->Update();

		typedef itk::ImageRegionIterator< ImageType2> IteratorType;	
		IteratorType it1 (reader->GetOutput(), reader->GetOutput()->GetRequestedRegion());
		
		for (it1.GoToBegin(); !it1.IsAtEnd(); ++it1){
			std::cout << it1.Get() << std::endl;
		}
		*/
	
		
		std::cout << "select the PCA to be computed" << std::endl;
	}

}
void MainWindow::recreateImageInstance(){
	
	//QStringList shapeFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (velocity fields)", "C:/0.Data/2. Statistical model/1. image/b. 130 images/5. validation/representation/original", "(*_velocity_field.mhd)");
	QStringList shapeFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (velocity fields)", "C:/0.Data/2. Statistical model/1. image/b. 130 images/5. validation/generalization/original", "(*_velocity_field.mhd)");
	QString shapeAverageFileName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/1. shape pca on 130 images/meanVelField.mhd");
	QStringList shapeEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors", "C:/0.Data/2. Statistical model/1. image/b. 130 images/1. shape pca on 130 images", "(*.mhd)");
		
	//QStringList intensityFileNames  = QFileDialog::getOpenFileNames(this, "Load intensities to recreate", "C:/0.Data/2. Statistical model/1. image/b. 130 images/5. validation/representation/original", "(*_rigid_warpedToRef.mhd)");
	QStringList intensityFileNames  = QFileDialog::getOpenFileNames(this, "Load intensities to recreate", "C:/0.Data/2. Statistical model/1. image/b. 130 images/5. validation/generalization/original", "(*_rigid_warpedToRef.mhd)"); 
	QString intensityAverageFileName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/2. intensity pca on 130 images/meanIntensity.mhd");
	QStringList intensityEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load intensity eigenvectors", "C:/0.Data/2. Statistical model/1. image/b. 130 images/2. intensity pca on 130 images", "(*.mhd)");
	
	QString combinedEVectorFileName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/3. combined pca on 130 images/combined model eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/3. combined pca on 130 images/combined model eigenvalues.txt");
	QString wFileName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/3. combined pca on 130 images/w and nOfImages.txt");
	//QString modesFileName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/5. validation/representation/every 10 modes/modes.txt");
	QString modesFileName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/5. validation/generalization/every 10 modes/modes.txt");
	//QString outputFolderName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/5. validation/representation/every 10 modes");
	QString outputFolderName = ("C:/0.Data/2. Statistical model/1. image/b. 130 images/5. validation/generalization/every 10 modes");
	

	double start = clock();

	PCAImages* PCA = new PCAImages;
	PCA->SetNumberOfModes(_modeInterval);
	
	PCA->SetInstanceShapeFileNames(shapeFileNames);
	PCA->SetShapeAverageFileName(shapeAverageFileName );
	PCA->SetShapeEVectorsFileNames(shapeEVectorFileNames);
			
	PCA->SetInstanceIntensityFileNames(intensityFileNames);
	PCA->SetIntensityAverageFileName(intensityAverageFileName);
	PCA->SetIntensityEVectorsFileNames(intensityEVectorFileNames);

	PCA->SetCombinedEVectorsFileName(combinedEVectorFileName);
	PCA->SetCombinedEValuesFileName(combinedEValuesFileName);
	
	PCA->SetWFileName(wFileName);
	PCA->SetModeNumbersFileName(modesFileName);
	PCA->SetOutputFolder(outputFolderName);
	PCA->InstanceRecreation();

	double end = clock();
	double total = (end-start)/CLOCKS_PER_SEC;
	std::cout << "Computation time for all instances: " << total << " sec. (" << total/60 << " min.)" << std::endl;
	
	/*
	// test data
	QStringList shapeFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (velocity fields)", "C:/0.Data/test data/pca validation/images/shape pca", "(*.mhd)");
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/images/shape pca/meanShape.mhd");
	QStringList shapeEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors", "C:/0.Data/test data/pca validation/images/shape pca", "(*.mhd)");
		
	QStringList intensityFileNames  = QFileDialog::getOpenFileNames(this, "Load intensities to recreate", "C:/0.Data/test data/pca validation/images/intensity pca", "(*.mhd)");
	QString intensityAverageFileName = ("C:/0.Data/test data/pca validation/images/intensity pca/meanIntensity.mhd");
	QStringList intensityEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load intensity eigenvectors", "C:/0.Data/test data/pca validation/images/intensity pca", "(*.mhd)");
	
	QString combinedEVectorFileName = ("C:/0.Data/test data/pca validation/images/combined pca/image combined pca eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/test data/pca validation/images/combined pca/image combined pca eigenvalues.txt");
	QString wFileName = ("C:/0.Data/test data/pca validation/images/combined pca/w and nOfImages.txt");
	QString modesFileName = ("C:/0.Data/test data/pca validation/images/modes.txt");

	QString outputFolderName = ("C:/0.Data/test data/pca validation/images");
	*/
	
	/*	
	// test data - 6 bones
	QStringList shapeFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (velocity fields)", "C:/0.Data/test data/pca validation/images/6 images/vfs", "(*.mhd)");
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/images/6 images/vfs/meanShape.mhd");
	QStringList shapeEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors", "C:/0.Data/test data/pca validation/images/6 images/vfs", "(*.mhd)");
		
	QStringList intensityFileNames  = QFileDialog::getOpenFileNames(this, "Load intensities to recreate", "C:/0.Data/test data/pca validation/images/6 images/warpedToReference", "(*.mhd)");
	QString intensityAverageFileName = ("C:/0.Data/test data/pca validation/images/6 images/warpedToReference/meanIntensity.mhd");
	QStringList intensityEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load intensity eigenvectors", "C:/0.Data/test data/pca validation/images/6 images/warpedToReference", "(*.mhd)");
	
	QString combinedEVectorFileName = ("C:/0.Data/test data/pca validation/images/6 images/combined/combined model eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/test data/pca validation/images/6 images/combined/combined model eigenvalues.txt");
	QString wFileName = ("C:/0.Data/test data/pca validation/images/6 images/combined/w and nOfImages.txt");
	QString modesFileName = ("C:/0.Data/test data/pca validation/images/6 images/recreated/modes.txt");

	QString outputFolderName = ("C:/0.Data/test data/pca validation/images/6 images/recreated");
	*/

	/*
	// only shape pca
	QStringList shapeFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (velocity fields)", "C:/0.Data/test data/pca validation/images/shape pca", "(*.mhd)");
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/images/shape pca/meanShape.mhd");
	QStringList shapeEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors", "C:/0.Data/test data/pca validation/images/shape pca", "(*.mhd)");
	QString outputFolderName = ("C:/0.Data/test data/pca validation/images/shape pca/recreated");
	*/
	
	/*
	QStringList shapeFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (velocity fields)", "C:/0.Data/test data/pca validation/images/6 images/vfs", "(*.mhd)");
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/images/6 images/vfs/meanShape.mhd");
	QStringList shapeEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors", "C:/0.Data/test data/pca validation/images/6 images/vfs", "(*.mhd)");
	QString outputFolderName = ("C:/0.Data/test data/pca validation/images/6 images/vfs/recreated");

	PCAImages* PCA = new PCAImages;
	PCA->SetInstanceShapeFileNames(shapeFileNames);
	PCA->SetShapeAverageFileName(shapeAverageFileName );
	PCA->SetShapeEVectorsFileNames(shapeEVectorFileNames);
	PCA->SetNumberOfModes(_modeInterval);
	PCA->SetOutputFolder(outputFolderName);
	PCA->RecreateShapePca();
	*/
	
	
}
void MainWindow::nOfImageInstances(int value){

	_nOfImageInstances = value;

}
void MainWindow::nOfImageModes(int value){

	_nOfImageModes = value;

}
void MainWindow::createImageInstances(){
	
	
	// true data
	QString shapeAverageFileName = QFileDialog::getOpenFileName(this, "Load average velocity field", "C:/0.Data/2011_FemaleMale_Study/3. shape pca/male", "(*.mhd)");
	QStringList shapeEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors", "C:/0.Data/2011_FemaleMale_Study/3. shape pca/male", "(*.mhd)");
	
	QString intensityAverageFileName = QFileDialog::getOpenFileName(this, "Load average intensity", "C:/0.Data/2011_FemaleMale_Study/4. intensity pca/male", "(*.mhd)");
	QStringList intensityEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load intensity eigenvectors", "C:/0.Data/2011_FemaleMale_Study/4. intensity pca/male", "(*.mhd)");
	
	QString combinedEVectorsFileName = QFileDialog::getOpenFileName(this, "Load combined eigenvectors", "C:/0.Data/2011_FemaleMale_Study/5. combined pca/male", "(*.txt)");
	QString combinedEValuesFileName = QFileDialog::getOpenFileName(this, "Load combined eigenvalues", "C:/0.Data/2011_FemaleMale_Study/5. combined pca/male", "(*.txt)");
	QString wFileName = QFileDialog::getOpenFileName(this, "Load w file name", "C:/0.Data/2011_FemaleMale_Study/5. combined pca/male", "(*.txt)");
	
	QString combinedWeightsFileName = QFileDialog::getOpenFileName(this, "Load combined weight", "C:/0.Data/2011_FemaleMale_Study/6. new instances/male", "(*.txt)");
	QString outputFolderName = QFileDialog::getExistingDirectory(this, "Load output folder", "C:/0.Data/2011_FemaleMale_Study/6. new instances/male");
	

	/*
	// test data - 6 images
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/images/6 images/vfs/meanShape.mhd");
	QStringList shapeEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors", "C:/0.Data/test data/pca validation/images/6 images/vfs", "(*.mhd)");
	
	QString intensityAverageFileName = ("C:/0.Data/test data/pca validation/images/6 images/warpedToReference/meanIntensity.mhd");
	QStringList intensityEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load intensity eigenvectors", "C:/0.Data/test data/pca validation/images/6 images/warpedToReference", "(*.mhd)");
	
	QString combinedEVectorsFileName = ("C:/0.Data/test data/pca validation/images/6 images/combined/combined model eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/test data/pca validation/images/6 images/combined/combined model eigenvalues.txt");
	QString wFileName = ("C:/0.Data/test data/pca validation/images/6 images/combined/w and nOfImages.txt");
	
	QString combinedWeightsFileName = ("C:/0.Data/test data/pca validation/images/6 images/newInstances/modes.txt");
	QString outputFolderName = ("C:/0.Data/test data/pca validation/images/6 images/newInstances");
	*/

	/*
	// test data
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/images/shape pca/meanShape.mhd");
	QStringList shapeEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load shape eigenvectors", "C:/0.Data/test data/pca validation/images/shape pca", "(*.mhd)");
	
	QString intensityAverageFileName = ("C:/0.Data/test data/pca validation/images/intensity pca/meanIntensity.mhd");
	QStringList intensityEVectorFileNames = QFileDialog::getOpenFileNames(this, "Load intensity eigenvectors", "C:/0.Data/test data/pca validation/images/intensity pca", "(*.mhd)");
	
	QString combinedEVectorsFileName = ("C:/0.Data/test data/pca validation/images/combined pca/image combined pca eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/test data/pca validation/images/combined pca/image combined pca eigenvalues.txt");
	QString wFileName = ("C:/0.Data/test data/pca validation/images/combined pca/w and nOfImages.txt");
	
	QString combinedWeightsFileName = ("C:/0.Data/test data/pca validation/images/createdInstances/modes.txt");
	QString outputFolderName = ("C:/0.Data/test data/pca validation/images/createdInstances");
	*/

	double start = clock();
	
	PCAImages* PCA = new PCAImages;
	PCA->SetNumberOfInstances(_nOfImageInstances);
	PCA->SetNumberOfModes(_nOfImageModes); 
	
	PCA->SetShapeAverageFileName(shapeAverageFileName );
	PCA->SetShapeEVectorsFileNames(shapeEVectorFileNames);

	PCA->SetIntensityAverageFileName(intensityAverageFileName);
	PCA->SetIntensityEVectorsFileNames(intensityEVectorFileNames);
			
	PCA->SetCombinedEVectorsFileName(combinedEVectorsFileName);
	PCA->SetCombinedEValuesFileName(combinedEValuesFileName);
	PCA->SetCombinedWeightsFileName(combinedWeightsFileName);

	PCA->SetWFileName(wFileName);
	PCA->SetOutputFolder(outputFolderName);
	PCA->InstanceCreation();

	double end = clock();
	double total = (end-start)/CLOCKS_PER_SEC;
	std::cout << "Computation time for all instances: " << total << " sec. (" << total/60 << " min.)" << std::endl;
	
}

/************************* MESH TOOL **************************/
// mesh creation
void MainWindow::loadImageMesh(){
	
	_fileNames = QFileDialog::getOpenFileNames(this, "Load files", " ", "Images (*.mhd);;Meshes (*.stl)");
}
void MainWindow::simplify(double value){

	if (value >= 0) 
		_simplifyValue = value;
	else 
		std::cout << "simplifyTargetReduction value invalid" << std::endl;

	//std::cout << "simplifyValue: " << _simplifyValue << std::endl;

}
void MainWindow::smooth(int value){

	if (value >= 0) 
		_smoothNofIteration = value;
	else 
		std::cout << "smoothNofIteration value invalid" << std::endl;

	//std::cout << "smoothNofIteration: " << _smoothNofIteration<< std::endl;
	
}
void MainWindow::createMesh(){

	
	// define if the files are images (.mhd) or meshes (.stl)
	bool fileType;

	if (_fileNames[0].endsWith(".mhd"))
		fileType = true;
	else if (_fileNames[0].endsWith(".stl"))
		fileType = false;
	else
		std::cout << "input files not supported" << std::endl;

	
	// for each file 
	for (int i=0; i<_fileNames.size(); i++){

		double start = clock();

		MeshSimplifyAndSmooth* simplifyAndSmooth = new MeshSimplifyAndSmooth;
	
		if (fileType == true){
			
			// load image
			ImageHandler* imageHandler = new ImageHandler;
			imageHandler->SetImageFileName(_fileNames[i]);
			std::cout << _fileNames[i].toAscii().data() << std::endl;
			imageHandler->MetafileReaderUpdate();
			imageHandler->SetThreshold(-800); // 
			imageHandler->MarchingCubesUpdate();
						
			// simplify and smooth
			simplifyAndSmooth->SetInput(imageHandler->GetMesh());
			simplifyAndSmooth->SetSimplifyValue(_simplifyValue);
			simplifyAndSmooth->SetSmoothIteration(_smoothNofIteration);
			simplifyAndSmooth->Update();
			
			// save
			_fileNames[i].replace(QString(".mhd"), QString(".stl"));
			std::cout << _fileNames[i].ascii() << std::endl;

			// cleaning up
			delete imageHandler;
		}
		 
		else if ( fileType == false){

			// load mesh
			vtkSTLReader* stlReader = vtkSTLReader::New();
			stlReader->SetFileName(_fileNames[i].ascii());
			std::cout << _fileNames[i].ascii() << std::endl;
			stlReader->Update();
			std::cout << "nodes: " << stlReader->GetOutput()->GetNumberOfPoints() << std::endl;
			std::cout << "elements: " <<  stlReader->GetOutput()->GetNumberOfCells() << std::endl;
			
			// simplify and smooth
			simplifyAndSmooth->SetInput(stlReader->GetOutput());
			simplifyAndSmooth->SetSimplifyValue(_simplifyValue);
			simplifyAndSmooth->SetSmoothIteration(_smoothNofIteration);
			simplifyAndSmooth->Update();

			// save
			_fileNames[i].replace(QString(".stl"), QString("_new.stl"));
			std::cout << _fileNames[i].ascii() << std::endl;
		
			stlReader->Delete();
		}

		// save
		std::cout << "mesh written" << std::endl;
		vtkSTLWriter* meshWriter = vtkSTLWriter::New();
		meshWriter->SetFileName(_fileNames[i].ascii());
		meshWriter->SetInput(simplifyAndSmooth->GetOutput());
		meshWriter->Update();
		
		// render
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(simplifyAndSmooth->GetOutput());
		renderingMesh->SetRenderer(_ren);
		//renderingMesh->Update();

		

		// cleaning up
		delete simplifyAndSmooth;
		meshWriter->Delete();
		delete renderingMesh;
		
		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "Computation time for mesh creation: " << total << " sec. about " << total/60 << " minutes" << std::endl;
		std::cout << std::endl;
	
	}
	

}
void MainWindow::loadSurfaceMesh(){
	
	_fileNames = QFileDialog::getOpenFileNames(this, "Load mesh files", " ", "Meshes (*.stl)");
}
void MainWindow::createVolumeMesh(){
	
	for (int i=0; i<_fileNames.size(); i++){

		double start = clock();
		
		std::cout << "run Netgen" << std::endl;		
		
		// NetGen
		QString NetGenPath = ("\"C:\\Program Files (x86)\\Netgen-4.9.13_x64\\bin\\netgen.exe\"");
		//QString NetGenPath = ("\"C:\\1.PhdProject\\tools\\appls\\netgen.exe\"");
		QString netgenCommand = (NetGenPath);
		
		netgenCommand.append(" -geofile=\""); // input
		netgenCommand.append(_fileNames[i]);
		netgenCommand.append("\"");
			
		netgenCommand.append (" -meshfile=\""); // output
		_fileNames[i].replace(QString(".stl"), QString(".inp"));
		netgenCommand.append(_fileNames[i]);

		netgenCommand.append ("\" -meshfiletype="); // output format
		netgenCommand.append("\"Abaqus Format\"");
		//netgenCommand.append(" -very fine"); // mesh granularity
		netgenCommand.append(" -batchmode");
		
		std::cout << netgenCommand.toAscii().data() << std::endl;
			
		QProcess myProcessNetgen;
		myProcessNetgen.start(netgenCommand.toAscii().data());
		myProcessNetgen.waitForFinished(200000);
		
		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "Computation time for NetGen: " << total << " sec" << std::endl;
	}

}
void MainWindow::meshJacobian(){
	
	_fileNames = QFileDialog::getOpenFileNames(this, "Load tetra mesh files", "C:/0.Data/1. Registered data/2a. morphed/volumes/abaqus inp", "Abaqus meshes (*.inp);;Ansys meshes (*.cdb)");

	std::cout << "calculate mesh element jacobian" << std::endl;
	
	double start = clock();

	// calculate jacobians
	MeshQuality* meshQuality = new MeshQuality;
	meshQuality->SetMeshFileNames(_fileNames);
	meshQuality->CalculateJacobian();

	double end = clock();
	double total = (end-start)/CLOCKS_PER_SEC;
	std::cout << "Computation time: " << total << " sec. about " << total/60 << " minutes" << std::endl;
	
	
}
void MainWindow::meshEdgeRatio(){
	
	_fileNames = QFileDialog::getOpenFileNames(this, "Load tetra mesh files", "C:/0.Data/1. Registered data/2a. morphed/volumes/abaqus inp", "Abaqus meshes (*.inp);;Ansys meshes (*.cdb)");

	std::cout << "calculate mesh element edge ratio" << std::endl;

	double start = clock();
	
	// calculate edge ratio
	MeshQuality* meshQuality = new MeshQuality;
	meshQuality->SetMeshFileNames(_fileNames);
	meshQuality->CalculateEdgeRatio();

	double end = clock();
	double total = (end-start)/CLOCKS_PER_SEC;
	std::cout << "Computation time: " << total << " sec. about " << total/60 << " minutes" << std::endl;
	
}
void MainWindow::meshAngle(){
	
	_fileNames = QFileDialog::getOpenFileNames(this, "Load tetra mesh files", "C:/0.Data/1. Registered data/2a. morphed/volumes/abaqus inp", "Abaqus meshes (*.inp);;Ansys meshes (*.cdb)");

	std::cout << "calculate mesh element min angle" << std::endl;

	double start = clock();
	
	// calculate jacobians
	MeshQuality* meshQuality = new MeshQuality;
	meshQuality->SetMeshFileNames(_fileNames);
	meshQuality->CalculateMinAngle();

	double end = clock();
	double total = (end-start)/CLOCKS_PER_SEC;
	std::cout << "Computation time: " << total << " sec. about " << total/60 << " minutes" << std::endl;
	
}

// morpher processing
void MainWindow::pickLandmarks(){

	// mesh
	_fileName = QFileDialog::getOpenFileName(this, "Load mesh", "C:/0.Data/1. Registered data/2. stl, txt and xml", "*.stl");
			
	std::cout << "mesh: " << _fileName.toAscii().data() << std::endl;
	vtkSTLReader* stlReader = vtkSTLReader::New();
	stlReader->SetFileName(_fileName);
	stlReader->Update();
	std::cout << "nodes: " << stlReader->GetOutput()->GetNumberOfPoints() << std::endl;
	std::cout << "elements: " <<  stlReader->GetOutput()->GetNumberOfCells() << std::endl;
	
	RenderingMesh* renderingMesh = new RenderingMesh;
	renderingMesh->SetMesh(stlReader->GetOutput());
	renderingMesh->SetRenderer(_ren);
	renderingMesh->Update();

	// picker
	std::cout << "pick landmarks (position the mouse on the point + press p)" << std::endl;
	_pointPicker->SetMesh(stlReader->GetOutput());
	_pointPicker->SetRenderer(_ren);
	double color[3];
	color[0] = 1.0; color[1] = 0.0; color[2] = 0.0;
	_pointPicker->SetColor(color);
	_pointPicker->SetPicker(_cellPicker);
	_coords = vtkPoints::New();
	_pointPicker->SetPoints(_coords);
	_pointPicker->Update();
	
}
void MainWindow::createXml(){
	
	QString referenceFileName = QFileDialog::getOpenFileName(this, "Load reference landmarks coordinates and id", "C:/0.Data/1. Registered data/2. stl, txt and xml", "*.txt");	
	
	_cellPicker->RemoveAllObservers();
	
	_coords = _pointPicker->GetPoints();
	for (int i=0; i<_coords->GetNumberOfPoints(); i++){
		double coord[3];
		_coords->GetPoint(i, coord);
	}

	// saving for visualization (to load in the renderMeshLandmarks)
	PointWriterXyz* pointWriter = new PointWriterXyz;
	pointWriter->SetInput(_coords);
	_fileName.replace(QString(".stl"), QString(".txt"));
	pointWriter->SetFileName(_fileName.toAscii().data());
	pointWriter->Update();
	std::cout << ".txt landmark coordinate file written" << std::endl; 
	
	// saving as .xml
	PointReaderXyzId* referenceLandmarkReader = new PointReaderXyzId;
	referenceLandmarkReader->SetFileName(referenceFileName);
	referenceLandmarkReader->Update();
	vtkDoubleArray* referenceLandmarksId = vtkDoubleArray::New();
	referenceLandmarksId = referenceLandmarkReader->GetIdVector();

	_fileName.replace(QString(".txt"), QString(".xml"));
	PointWriterMorpherXml* xmlWriter = new PointWriterMorpherXml;
	xmlWriter->SetFileName(_fileName);
	xmlWriter->SetInput(_coords);
	xmlWriter->SetIdVector(referenceLandmarksId);
	referenceFileName.replace(QString(".txt"), QString(".cdb"));
	xmlWriter->SetReferenceFileName(referenceFileName);
	xmlWriter->Update();
	std::cout << ".xml file written" << std::endl; 

	// cleaning up
	delete pointWriter;
	delete referenceLandmarkReader;
	delete xmlWriter;
	_coords->Delete();
	
}
void MainWindow::createXmlFromFiles(){

	_fileNames = QFileDialog::getOpenFileNames(this, "Load landmarks coordinates", " ", "Landmarks (*.txt)");
	
	QString referenceFileName = QFileDialog::getOpenFileName(this, "Load reference landmarks coordinates and id", " ", "Landmarks and Id (*.txt)");	
	
	PointWriterMorpherXml* xmlWriter = new PointWriterMorpherXml;
	
	// reference IDs
	PointReaderXyzId* referenceLandmarkReader = new PointReaderXyzId;
	referenceLandmarkReader->SetFileName(referenceFileName);
	referenceLandmarkReader->Update();
	vtkDoubleArray* referenceLandmarksId = vtkDoubleArray::New();
	referenceLandmarksId = referenceLandmarkReader->GetIdVector();
	xmlWriter->SetIdVector(referenceLandmarksId);
	referenceFileName.replace(QString(".txt"), QString(".cdb"));
	xmlWriter->SetReferenceFileName(referenceFileName);
	
	PointReaderXyz* landmarksReader = new PointReaderXyz;	
	
	for (int i=0; i<_fileNames.size(); i++){
		
		// landmarks coordinates
		std::cout << _fileNames[i].toAscii().data() << std::endl; 
		landmarksReader->SetFileName(_fileNames[i]);
		landmarksReader->Update();
		vtkPoints* landmarks = vtkPoints::New();
		landmarks = landmarksReader->GetOutput();
		xmlWriter->SetInput(landmarks);
	
		// file creation
		_fileNames[i].replace(QString(".txt"), QString(".xml"));
		xmlWriter->SetFileName(_fileNames[i]);
		xmlWriter->Update();
	}

	// cleaning up
	delete xmlWriter;
	delete referenceLandmarkReader;
	delete landmarksReader;

	/*
	// reference
	QString referenceFileName = QFileDialog::getOpenFileName(this, "Load reference landmarks coordinates and id", " ", "Landmarks and Id (*.txt)");	
	PointReaderXyzId* referenceLandmarkReader = new PointReaderXyzId;
	referenceLandmarkReader->SetFileName(referenceFileName);
	referenceLandmarkReader->Update();
	vtkIntArray* referenceLandmarksId = vtkIntArray::New();
	referenceLandmarksId = referenceLandmarkReader->GetIdVector();
	
	//xml writer
	PointWriterMorpherXml* xmlWriter = new PointWriterMorpherXml;
	xmlWriter->SetIdVector(referenceLandmarksId);
	referenceFileName.replace(QString(".txt"), QString(".cdb"));
	xmlWriter->SetReferenceFileName(referenceFileName);

	// landmarks
	_fileNames = QFileDialog::getOpenFileNames(this, "Load landmarks coordinates", " ", "Landmarks (*.txt)");
	
	PointReaderXyz* landmarksReader = new PointReaderXyz;	
	vtkPoints* landmarks = vtkPoints::New();
	for (int i=0; i<_fileNames.size(); i++){
		
		// landmarks coordinates
		std::cout << _fileNames[i].toAscii().data() << std::endl; 
		landmarksReader->SetFileName(_fileNames[i]);
		landmarksReader->Update();
		landmarks = landmarksReader->GetOutput();
		landmarks->GetNumberOfPoints();
		
		vtkPoints* selectedLandmarks = vtkPoints::New();
		for (int i=0; i<10; i++){
			double pt[3];
			landmarks->GetPoint(landmarks->GetNumberOfPoints()-10+i, pt);
			selectedLandmarks->InsertNextPoint(pt);
		}

		xmlWriter->SetInput(selectedLandmarks);
	
		// txt overwriting
		PointWriterXyz* pointWriter = new PointWriterXyz;
		pointWriter->SetFileName(_fileNames[i]);
		pointWriter->SetInput(selectedLandmarks);
		pointWriter->Update();

		// xml writing
		_fileNames[i].replace(QString(".txt"), QString(".xml"));
		xmlWriter->SetFileName(_fileNames[i]);
		xmlWriter->Update();
		std::cout << "files written" << std::endl;
		
		

	}

	// cleaning up
	delete xmlWriter;
	delete referenceLandmarkReader;
	delete landmarksReader;	
	*/
}
void MainWindow::runMorpher(){

	_fileNames = QFileDialog::getOpenFileNames(this, "Load .xml files", " ", "xml files (*.xml)");
		
	QString AnsysMorpherPath = ("\"C:\\1.PhdProject\\tools\\appls\\ansys morpher\\exe\\Morpher.exe\"");

	double start = clock();

	for (int i=0; i<_fileNames.size(); i++){

		// run ansys morpher
		std::cout << "mesh to morph: " << _fileNames[i].toAscii().data() << std::endl;
		std::cout << "running ansys mesh morpher" << std::endl;
		
		QStringList arguments;
		arguments << _fileNames[i].toAscii();
				
		QProcess myProcess;
		myProcess.start(AnsysMorpherPath, arguments);
		myProcess.waitForFinished(2000000);
		std::cout << "mesh morphed" << std::endl;
	}

	double end = clock();
	double total = (end-start)/CLOCKS_PER_SEC;
	std::cout << "computation time of morphing " << _fileNames.size() << " meshes: " << total << " sec. about " << total/60 << " minutes" << std::endl;
}
void MainWindow::surfaceValidation(){

	QStringList volumeMeshFileNames = QFileDialog::getOpenFileNames(this, "Load registered meshes", "C:/0.Data/1. Registered data/2a. morphed/volumes/ansys cdb", "*.cdb");
	QStringList stlMeshFileNames = QFileDialog::getOpenFileNames(this, "Load original meshes", "C:/0.Data/1. Registered data/2. stl, txt and xml", "*.stl");

	RegistrationValidationSurfaceMesh* validation = new RegistrationValidationSurfaceMesh;
	validation->SetFlag(0);
	validation->SetRegisteredFileNames(volumeMeshFileNames);
	validation->SetOriginalFileNames(stlMeshFileNames);
	validation->Update();
		
	delete validation;
}
void MainWindow::volumeValidation(){

	// reference volume .cdb mesh file name
	QString referenceMeshFileName = QFileDialog::getOpenFileName(this, "Load reference volume mesh", "C:/0.Data/1. Registered data/0. ref 00705", "*.cdb");
	// reference mask file name
	QString referenceMaskFileName = QFileDialog::getOpenFileName(this, "Load reference masked image", "C:/0.Data/1. Registered data/0. ref 00705", "*Classified-SLICES.mhd");
	// moving volume .cdb mesh file names
	QStringList movingMeshFileNames = QFileDialog::getOpenFileNames(this, "Load Ansys Morpher volume meshes", "C:/0.Data/1. Registered data/2a. morphed/volumes/ansys cdb/135 correspondent to mask - registration validation", "*.cdb");
	// moving mask file names
	QStringList movingMaskFileNames = QFileDialog::getOpenFileNames(this, "Load the masked images correspondent to the meshes", "C:/0.Data/0. Original data/3. femur left masks/0a. femur left masks - three layers", "*.mhd");

	// volume validation
	RegistrationValidationVolumeMesh* validation = new RegistrationValidationVolumeMesh;
	validation->SetReferenceMeshFileName(referenceMeshFileName);
	validation->SetReferenceMaskFileName(referenceMaskFileName);
	validation->SetMovingMeshFileNames(movingMeshFileNames);
	validation->SetMovingMaskFileNames(movingMaskFileNames);
	validation->Update();

	delete validation;
	
	
}

// statistical model
void MainWindow::extractCoordinatesAndIntensities(){
	
	// when reading both mesh and image (nodes and intensities extraction)
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load meshes","C:/0.Data/1. Registered data/2a. morphed/volumes/ansys cdb", "Ansys quadratic volume meshes (*.cdb);;Abaqus quadratic volume meshes (*.inp)");
	QStringList imageFileNames = QFileDialog::getOpenFileNames(this, "Load images","C:/0.Data/0. Original data/1. femur left mhd", "Images (*.mhd)");
		
	for (int i=0; i<meshFileNames.size(); i++){
		
		vtkPolyData* mesh = vtkPolyData::New();

		// ansys mesh
		if (meshFileNames[i].endsWith(".cdb")){ // Ansys quadratic volume meshes
			MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
			meshReader->SetFileName(meshFileNames[i]);
			std::cout << meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			mesh->DeepCopy(meshReader->GetOutput());
			delete meshReader;		
		}	
		
		// abaqus mesh
		if (meshFileNames[i].endsWith(".inp")){
			MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
			meshReader->SetFileName(meshFileNames[i]);
			std::cout << meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			mesh->DeepCopy(meshReader->GetOutput());
			delete meshReader;
		}
		
		// read correspondent image (for stryker data strings)
		QString temp = meshFileNames[i];
		if (temp.lastIndexOf("/") == -1){
			temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
			temp.remove(0, temp.lastIndexOf("\\")+1);
		}
		else {
			temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
			temp.remove(0, temp.lastIndexOf("/")+1);
		}
		ImageHandlerFloat* imageHandler = new ImageHandlerFloat;
		for (int y=0; y<imageFileNames.size(); y++){
			if (imageFileNames[y].contains(temp)){
				std::cout << "image: " << imageFileNames[y].ascii() << std::endl;
				imageHandler->SetImageFileName(imageFileNames[y].ascii());
				imageHandler->MetafileReaderUpdate();
			}
		}
		// image extrusion
		std::cout << "extruding and getting grey levels" << std::endl;
		imageHandler->Extrusion();

		// extract nodes grey levels
		FemAssignerNodes* greyLevelsExtractor = new FemAssignerNodes;
		greyLevelsExtractor->SetMesh(mesh);
		greyLevelsExtractor->SetImage(imageHandler->GetImage()); // to uncomment after changing the image type in FemAssigner.h
		greyLevelsExtractor->GreyLevelAssignmentUpdate();
		vtkDoubleArray* greyLevels = vtkDoubleArray::New();
		greyLevels = greyLevelsExtractor->GetGreyLevels();

		std::cout << "writing files" << std::endl;
		// nodes filename
		temp = meshFileNames[i];
		if (temp.lastIndexOf("/") == -1){
			temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
			temp.append("_nodes.txt");
		}
		else {
			temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
			temp.append("_nodes.txt");
		}
		std::cout << temp.toAscii().data() << std::endl;
		// nodes to vnl_vector
		vnl_vector<double> nodes;
		nodes.set_size(mesh->GetNumberOfPoints()*3);
		for (int i=0; i<mesh->GetNumberOfPoints(); i++){
			double pt[3];
			mesh->GetPoint(i,pt);
			nodes(i*3) = pt[0];
			nodes(i*3+1) = pt[1];
			nodes(i*3+2) = pt[2];
		}
		// write vnl
		VnlWriterVector* vectorWriter = new VnlWriterVector;
		vectorWriter->SetFileName(temp);
		vectorWriter->SetVnlVector(nodes);
		vectorWriter->Update();

		// intensity file name
		temp.replace(QString("_nodes.txt"), QString("_intensities.txt"));
		std::cout << temp.toAscii().data() << std::endl;
		// intensity to vnl_vector
		vnl_vector<double> intensities;
		intensities.set_size(greyLevels->GetNumberOfTuples());
		for (int i=0; i<greyLevels->GetNumberOfTuples(); i++){
			intensities(i) = greyLevels->GetValue(i) * 1000; // grey level are divided by 1000 in the FemAssigner class
		}
		// write vnl
		vectorWriter->SetFileName(temp);
		vectorWriter->SetVnlVector(intensities);
		vectorWriter->Update();

		// cleaning
		delete imageHandler;
		delete greyLevelsExtractor;
		delete vectorWriter;

		std::cout << std::endl;
	}
	std::cout << "files written" << std::endl;
	
		
	/*
	// when reading only mesh (node extraction)
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load aligned meshes","C:/0.Data/1. Registered data/2a. morphed/aligned volumes", "Ansys quadratic volume meshes (*.cdb)");
		
	for (int i=0; i<meshFileNames.size(); i++){
			
			// read mesh (ansys .cdb)
			MeshReaderMorpherVolume* reader = new MeshReaderMorpherVolume;
			std::cout << "mesh: " << meshFileNames[i].ascii() << std::endl;
			reader->SetFileName(meshFileNames[i]);
			reader->Update();
			
			// nodes filename
			std::cout << "writing file" << std::endl;
			QString temp = meshFileNames[i];
			temp = meshFileNames[i];
			if (temp.lastIndexOf("/") == -1){
				temp.remove(temp.lastIndexOf("\\")+6,temp.size()-1);
				temp.append("_nodes_aligned.txt");
			}
			else {
				temp.remove(temp.lastIndexOf("/")+6,temp.size()-1);
				temp.append("_nodes_aligned.txt");
			}
			std::cout << temp.toAscii().data() << std::endl;
			// nodes to vnl_vector
			vnl_vector<double> nodes;
			nodes.set_size(reader->GetOutput()->GetNumberOfPoints()*3);
			for (int i=0; i<reader->GetOutput()->GetNumberOfPoints(); i++){
				double pt[3];
				reader->GetOutput()->GetPoint(i,pt);
				nodes(i*3) = pt[0];
				nodes(i*3+1) = pt[1];
				nodes(i*3+2) = pt[2];
			}
			// write vnl
			VnlWriterVector* vectorWriter = new VnlWriterVector;
			vectorWriter->SetFileName(temp);
			vectorWriter->SetVnlVector(nodes);
			vectorWriter->Update();

			// cleaning
			delete reader;
			delete vectorWriter;

			std::cout << std::endl;
	}
	std::cout << "files written" << std::endl;
	*/
	
}
void MainWindow::extractNodesAndElements(){

	// when reading only mesh (node extraction)
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load meshes","C:/0.Data/2011_FemaleMale_Study/6. new instances/female/meshes/inp", "Abaqus quadratic volume mesh (*.inp);;Ansys quadratic volume meshes (*.cdb)");
		
	for (int i=0; i<meshFileNames.size(); i++){
			
		vtkPolyData* mesh = vtkPolyData::New();
		
		// abaqus mesh
		if (meshFileNames[i].endsWith(".inp")){
			MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
			meshReader->SetFileName(meshFileNames[i]);
			std::cout << meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			mesh->DeepCopy(meshReader->GetOutput());
			delete meshReader;
		}
		
		// ansys mesh
		if (meshFileNames[i].endsWith(".cdb")){ // Ansys quadratic volume meshes
			MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
			meshReader->SetFileName(meshFileNames[i]);
			std::cout << meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			mesh->DeepCopy(meshReader->GetOutput());
			delete meshReader;		
		}
	
		// write nodes 
		QString temp = meshFileNames[i];
		temp.remove(temp.lastIndexOf("."), temp.size()-1);
		temp.append("_nodes.txt");
		std::cout << temp.toAscii().data() << std::endl;
		PointWriterXyz* pointWriter = new PointWriterXyz;
		pointWriter->SetInput(mesh->GetPoints());
		pointWriter->SetFileName(temp);
		pointWriter->Update();
		delete pointWriter;
		
		// write elements
		temp.replace(QString("_nodes.txt"), QString("_elements.txt"));
		std::cout << temp.toAscii().data() << std::endl;
		
		QFile outFile(temp.toAscii().data());
		outFile.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream file(&outFile);
		vtkGenericCell* cell = vtkGenericCell::New();
		vtkIdList* idList;
		for (int i=0; i<mesh->GetNumberOfCells(); i++){
			mesh->GetCell(i, cell);
			idList = cell->GetPointIds();
			for (int a=0; a<10; a++){ // 10 nodes
				if (a == 9)
					file << idList->GetId(a)+1 << endl;
				else
				file << idList->GetId(a)+1 << " " ;
				
			}
		}
		outFile.close();

		/*
		// write first node (for rendering)
		vtkPoints* firstNode = vtkPoints::New();
		for (int i=0; i<mesh->GetNumberOfCells(); i++){
			// get the element
			mesh->GetCell(i, cell);
			idList = cell->GetPointIds();
			// get the first node
			double pt[3];
			mesh->GetPoint(idList->GetId(0), pt);
			firstNode->InsertNextPoint(pt);
		}

		PointWriterXyz* writer = new PointWriterXyz;
		temp.replace(QString("_elements.txt"), QString("_firstNodeOfEachElement.txt"));
		std::cout << temp.toAscii().data() << std::endl;
		writer->SetFileName(temp);
		writer->SetInput(firstNode);
		writer->Update();
		delete writer;
		*/

	}
	
	std::cout << "files written" << std::endl;
	

}
void MainWindow::meshAlignment(){

	/*
	_fileNames = QFileDialog::getOpenFileNames(this, "Load meshes"," ","Ansys quadratic volume meshes (*.cdb)");

	RegistrationProcrustesAlignment* procrustes = new RegistrationProcrustesAlignment;
	procrustes->SetFileNames(_fileNames);
	procrustes->SetRegistrationFlag(0);
	procrustes->Update();
	*/

	
	///////////////////////////////////////////
	/////////////////// ICP ///////////////////
	///////////////////////////////////////////
	QString sourceFileName = QFileDialog::getOpenFileName(this, "Load source mesh"," ","Abaqus quadratic volume meshes (*.inp)");
	//QString targetFileName = QFileDialog::getOpenFileName(this, "Load target mesh"," ","Ansys quadratic volume meshes (*.cdb)");
	QString targetFileName = QFileDialog::getOpenFileName(this, "Load target mesh"," ","Abaqus quadratic volume meshes (*.inp)");
	std::cout << "source: " << sourceFileName.toAscii().data() << std::endl;
	std::cout << "target: " << targetFileName.toAscii().data() << std::endl;



	// read source
	MeshReaderAbaqus* readerSource = new MeshReaderAbaqus;
	readerSource->SetFileName(sourceFileName);
	readerSource->Update();

	// read target
	MeshReaderAbaqus* readerTarget = new MeshReaderAbaqus;
	//MeshReaderMorpherVolume* readerTarget = new MeshReaderMorpherVolume;
    readerTarget->SetFileName(targetFileName);
	readerTarget->Update();
		
	// extract the outer surface for alignment
	//MeshExtractOuterSurface* sourceExtractor = new MeshExtractOuterSurface;
	//sourceExtractor->SetVolumeMesh(readerSource->GetOutput());
	//sourceExtractor->Update();

	//MeshExtractOuterSurface* targetExtractor = new MeshExtractOuterSurface;
	//targetExtractor->SetVolumeMesh(readerTarget->GetOutput());
	//targetExtractor->Update();

	// align
	RegistrationProcrustesAlignment* ICP = new RegistrationProcrustesAlignment;
	ICP->SetSourceMesh(readerSource->GetOutput());
	ICP->SetTargetMesh(readerTarget->GetOutput());
	//ICP->SetSourceMesh(sourceExtractor->GetSurfaceMesh());
	//ICP->SetTargetMesh(targetExtractor->GetSurfaceMesh());
	ICP->IterativeClosestPoints();

	// save aligned
	MeshWriterAbaqus* writer = new MeshWriterAbaqus;
	sourceFileName.replace(QString(".inp"), QString("_aligned.inp"));
	writer->SetFileName(sourceFileName);
	writer->MeshOn();
	writer->SetMesh(ICP->GetAlignedMesh());
	writer->Update();
	std::cout << ".inp file written" << std::endl;


	// visualization aligned
	MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
	surfaceExtractor->SetVolumeMesh(ICP->GetAlignedMesh());
	surfaceExtractor->Update();
	double color[3];
	color[0] = 232.0/255.0; color[1] = 184.0/255.0; color[2] = 45.0/255.0; // bone color
	RenderingMesh* renderingMesh = new RenderingMesh;
	renderingMesh->SetColor(color);
	renderingMesh->SetMesh(surfaceExtractor->GetSurfaceMesh());
	renderingMesh->SetRenderer(_ren);
	renderingMesh->Update();

	// visualization target
	MeshExtractOuterSurface* surfaceExtractorTarget = new MeshExtractOuterSurface;
	surfaceExtractorTarget->SetVolumeMesh(readerTarget->GetOutput());
	surfaceExtractorTarget->Update();
	color[0] = 0.0/255.0; color[1] = 255.0/255.0; color[2] = 0.0/255.0; // bone color
	RenderingMesh* renderingMeshTarget = new RenderingMesh;
	renderingMeshTarget->SetColor(color);
	renderingMeshTarget->SetMesh(surfaceExtractorTarget->GetSurfaceMesh());
	renderingMeshTarget->SetRenderer(_ren);
	renderingMeshTarget->Update();
	
	delete readerSource;
	delete readerTarget;
	delete ICP;
	delete writer;
	delete surfaceExtractor;
	delete renderingMesh;
	delete surfaceExtractorTarget;
	delete renderingMeshTarget;

}
void MainWindow::shapeMeshPCA(){

	_PCAflag = 1;
	
}
void MainWindow::instensityMeshPCA(){

	_PCAflag = 2;
	
}
void MainWindow::combinedMeshPCA(){

	_PCAflag = 3;
	

}
void MainWindow::calculateMeshPCA(){

	// shape PCA
	if (_PCAflag == 1){

		std::cout << std::endl;
		std::cout << "compute shape PCA" << std::endl;
		std::cout << "be sure that you are loading aligned meshes with correspondent nodes!" << std::endl;
		
		QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load meshes","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/0. meshes/aligned","Ansys quadratic volume meshes (*.cdb);;Abaqus quadratic volume meshes (*.inp);;Surface linear meshes (*.stl)");
		
		double start = clock();
		
		PCAMesh* PCA = new PCAMesh;
		PCA->SetMeshFileNames(fileNames);
		PCA->ShapePCA();

		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "Computation time: " << total << " sec. (" << total/60 << " min.)" << std::endl;
	}

	// instensity PCA
	else if (_PCAflag == 2){
		
		std::cout << std::endl;
		std::cout << "compute intensity PCA" << std::endl;
		std::cout << "be sure that you are loading meshes and images that are correspondent!" << std::endl;
		
		QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load meshes","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/0. meshes/original", "Ansys quadratic volume meshes (*.cdb);;Abaqus quadratic volume meshes (*.inp);;Surface linear meshes (*.stl)");
		QStringList imageFileNames = QFileDialog::getOpenFileNames(this, "Load images","C:/0.Data/0. Original data/1. femur left mhd", "Images (*.mhd)");
		
		double start = clock();
		
		PCAMesh* PCA = new PCAMesh;
		PCA->SetMeshFileNames(meshFileNames);
		PCA->SetImageFileNames(imageFileNames);
		PCA->IntensityPCA();
		
		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "Computation time: " << total << " sec. (" << total/60 << " min.)" << std::endl;

	}

	// combined PCA
	else if (_PCAflag == 3){
		
		std::cout << std::endl;
		std::cout << "compute combined PCA" << std::endl;
		
		QString shapeDatasetfileName = QFileDialog::getOpenFileName(this, "Load shape dataset file","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/1. shape pca", "*.txt");
		QString shapeAveragefileName = QFileDialog::getOpenFileName(this, "Load shape average file","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/1. shape pca", "*.txt");
		QString shapeEvaluesfileName = QFileDialog::getOpenFileName(this, "Load shape eigenvalue file","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/1. shape pca", "*.txt");
		QString shapeEvectorsfileName = QFileDialog::getOpenFileName(this, "Load shape eigenvector file","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/1. shape pca", "*.txt");
		
		QString intensityDatasetfileName = QFileDialog::getOpenFileName(this, "Load intensity dataset file","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/2. intensity pca", "*.txt");
		QString intensityAveragefileName = QFileDialog::getOpenFileName(this, "Load intensity average file","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/2. intensity pca", "*.txt");
		QString intensityEvaluesfileName = QFileDialog::getOpenFileName(this, "Load intensity eigenvalues file","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/2. intensity pca", "*.txt");
		QString intensityEvectorsfileName = QFileDialog::getOpenFileName(this, "Load intensity eigenvector file","C:/0.Data/2011_FemaleMale_Study/0. meshes/male/2. intensity pca", "*.txt"); 
		

		double start = clock();

		PCAMesh* PCA = new PCAMesh;

		PCA->SetShapeDatasetCoordinatesFileName(shapeDatasetfileName);
		PCA->SetShapeAverageFileName(shapeAveragefileName);
		PCA->SetShapeEValuesFileName(shapeEvaluesfileName);
		PCA->SetShapeEVectorsFileName(shapeEvectorsfileName);

		PCA->SetIntensityDatasetGreyLevelFileName(intensityDatasetfileName);
		PCA->SetIntensityAverageFileName(intensityAveragefileName);
		PCA->SetIntensityEValuesFileName(intensityEvaluesfileName);
		PCA->SetIntensityEVectorsFileName(intensityEvectorsfileName);
		
		PCA->CombinedPCA();

		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "Computation time: " << total << " sec. (" << total/60 << " min.)" << std::endl;
	}

	else if (_PCAflag != 1 && _PCAflag != 2 && _PCAflag != 3){

		/*
		// testMatrixPCA
		QString fileName = QFileDialog::getOpenFileName(this, "Load file","C:/0.Data/test data/pca validation/images/matlab/mesh code", "*.txt");
		PCAMesh* PCA2 = new PCAMesh;
		PCA2->SetShapeDatasetCoordinatesFileName(fileName);
		PCA2->TestMatrixPCA();
		std::cout << "done" << std::endl;
		*/
		std::cout << "select the PCA to be computed" << std::endl;
			
	}

}
void MainWindow::recreateMeshInstance(){
	

	//QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (aligned meshes)", "C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/5. validation/generalization/original", "(*.cdb)");
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (aligned meshes)", "C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/5. validation/representation/original", "(*.cdb)");
	QString shapeAverageFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/1. shape pca on 130 meshes/mesh shape pca average.cdb");
	QString shapeEVectorFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/1. shape pca on 130 meshes/mesh shape pca eigenvectors.txt");
	
	//QStringList intensityFileNames  = QFileDialog::getOpenFileNames(this, "Load intensities to recreate", "C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/5. validation/generalization/original", "(*.txt)");
	QStringList intensityFileNames  = QFileDialog::getOpenFileNames(this, "Load intensities to recreate", "C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/5. validation/representation/original", "(*.txt)");
	QString intensityAverageFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/2. intensity pca on 130 meshes/mesh intensity pca average.txt");
	QString intensityEVectorFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/2. intensity pca on 130 meshes/mesh intensity pca eigenvectors.txt");
	
	QString combinedEVectorFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/3. combined pca on 130 meshes/mesh combined pca eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/3. combined pca on 130 meshes/mesh combined pca eigenvalues.txt");
	QString wFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/3. combined pca on 130 meshes/w and nOfDataset.txt");
	//QString modesFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/5. validation/generalization/every 10 modes/modes.txt");
	QString modesFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/5. validation/representation/modes.txt");

	//QString outputFolderName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/5. validation/generalization");
	QString outputFolderName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/5. validation/representation");
	
	
	PCAMesh* PCA = new PCAMesh;
		
	PCA->SetInstanceShapeFileNames(meshFileNames);
	PCA->SetShapeAverageFileName(shapeAverageFileName );
	PCA->SetShapeEVectorsFileName(shapeEVectorFileName);
			
	PCA->SetInstanceIntensityFileNames(intensityFileNames);
	PCA->SetIntensityAverageFileName(intensityAverageFileName);
	PCA->SetIntensityEVectorsFileName(intensityEVectorFileName);

	PCA->SetCombinedEVectorsFileName(combinedEVectorFileName);
	PCA->SetCombinedEValuesFileName(combinedEValuesFileName);
	
	PCA->SetWFileName(wFileName);
	PCA->SetModeNumbersFileName(modesFileName);
	PCA->SetOutputFolder(outputFolderName);
	PCA->InstanceRecreation();
	
/*
	// test with syntetic data
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate", "C:/0.Data/test data/pca validation/mesh/test data");
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/mesh/test data/shape average.txt");
	QString shapeEVectorFileName = ("C:/0.Data/test data/pca validation/mesh/test data/shape eigenvectors.txt");
	
	QStringList intensityFileNames = QFileDialog::getOpenFileNames(this, "Load intensity to recreate", "C:/0.Data/test data/pca validation/mesh/test data");
	QString intensityAverageFileName = ("C:/0.Data/test data/pca validation/mesh/test data/intensity average.txt");
	QString intensityEVectorFileName = ("C:/0.Data/test data/pca validation/mesh/test data/intensity eigenvectors.txt");
	
	QString combinedEVectorFileName = ("C:/0.Data/test data/pca validation/mesh/test data/combined eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/test data/pca validation/mesh/test data/combined eigenvalues.txt");

	QString wFileName = ("C:/0.Data/test data/pca validation/mesh/test data/w and nOfMesh.txt");
	QString modesFileName = ("C:/0.Data/test data/pca validation/mesh/true data/modes.txt");
	QString outputFolderName = ("C:/0.Data/test data/pca validation/mesh/test data");
	*/
	
	
	// test with 6 meshes
	/*
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (aligned meshes)", "C:/0.Data/test data/pca validation/mesh/true data/1. shape dataset", "(*.cdb)");
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/mesh/true data/mesh shape pca average.cdb");
	QString shapeEVectorFileName = ("C:/0.Data/test data/pca validation/mesh/true data/mesh shape pca eigenvectors.txt");
	
	QStringList intensityFileNames = QFileDialog::getOpenFileNames(this, "Load intensity to recreate", "C:/0.Data/test data/pca validation/mesh/true data/2. intensity dataset", "(*.txt)");
	QString intensityAverageFileName = ("C:/0.Data/test data/pca validation/mesh/true data/mesh intensity pca average.txt");
	QString intensityEVectorFileName = ("C:/0.Data/test data/pca validation/mesh/true data/mesh intensity pca eigenvectors.txt");
	
	QString bFileName = ("C:/0.Data/test data/pca validation/mesh/true data/b.txt");
	QString combinedEVectorFileName = ("C:/0.Data/test data/pca validation/mesh/true data/mesh combined pca eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/test data/pca validation/mesh/true data/mesh combined pca eigenvalues.txt");
		
	QString wFileName = ("C:/0.Data/test data/pca validation/mesh/true data/w and nOfMesh.txt");
	QString modesFileName = ("C:/0.Data/test data/pca validation/mesh/true data/modes.txt");
	
	QString outputFolderName = ("C:/0.Data/test data/pca validation/mesh/true data");
	*/

	/*
	// test for comparing with images
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate", "C:/0.Data/test data/pca validation/images/matlab/mesh code");
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code/shape average.txt");
	QString shapeEVectorFileName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code/shape eigenvectors.txt");
	
	QStringList intensityFileNames = QFileDialog::getOpenFileNames(this, "Load intensity to recreate", "C:/0.Data/test data/pca validation/images/matlab/mesh code");
	QString intensityAverageFileName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code/intensity average.txt");
	QString intensityEVectorFileName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code/intensity eigenvectors.txt");
	
	QString combinedEVectorFileName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code/combined eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code/combined eigenvalues.txt");

	QString wFileName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code/w and nOfMesh.txt");
	QString modesFileName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code/modes.txt");
	QString outputFolderName = ("C:/0.Data/test data/pca validation/images/matlab/mesh code");
	*/
	/*
	// test only shape
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load shapes to recreate (aligned meshes)", "C:/0.Data/test data/pca validation/mesh/true data/1. shape dataset", "(*.cdb)");
	QString shapeAverageFileName = ("C:/0.Data/test data/pca validation/mesh/true data/mesh shape pca average.cdb");
	QString shapeEVectorFileName = ("C:/0.Data/test data/pca validation/mesh/true data/mesh shape pca eigenvectors.txt");
	
	QString outputFolderName = ("C:/0.Data/test data/pca validation/mesh/true data");

	PCAMesh* PCA = new PCAMesh;
	PCA->SetNumberOfModes(_modeInterval);
	PCA->SetInstanceShapeFileNames(meshFileNames);
	PCA->SetShapeAverageFileName(shapeAverageFileName );
	PCA->SetShapeEVectorsFileName(shapeEVectorFileName);
	PCA->SetOutputFolder(outputFolderName);
	PCA->ShapeInstanceRecreation();
	*/
	
}
void MainWindow::nOfMeshInstances(int value){

	_nOfMeshInstances = value;
}
void MainWindow::nOfMeshModes(int value){

	_nOfMeshModes = value;

}
void MainWindow::createMeshInstances(){

	/*
	QString shapeAverageFileName = ("C:/0.Data/test data/debug image-mesh pca/meshes/shape average.txt");
	QString shapeEVectorFileName = ("C:/0.Data/test data/debug image-mesh pca/meshes/image shape eigenvectors.txt");
	
	QString intensityAverageFileName = ("C:/0.Data/test data/debug image-mesh pca/meshes/intensity average.txt");
	QString intensityEVectorFileName = ("C:/0.Data/test data/debug image-mesh pca/meshes/image intensity eigenvectors.txt");
	
	QString combinedEVectorFileName = ("C:/0.Data/test data/debug image-mesh pca/meshes/image combined pca eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/test data/debug image-mesh pca/meshes/image combined pca eigenvalues.txt");
	QString combinedWeightsFileName = ("C:/0.Data/test data/debug image-mesh pca/meshes/weights.txt");
	
	QString wFileName = ("C:/0.Data/test data/debug image-mesh pca/images/w and nOfImages.txt");
	
	QString outputFolderName = ("C:/0.Data/test data/debug image-mesh pca/meshes/to use mesh code with image values");
	*/

	
	QString shapeAverageFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/1. shape pca on 130 meshes/mesh shape pca average.cdb");
	QString shapeEVectorFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/1. shape pca on 130 meshes/mesh shape pca eigenvectors.txt");
	
	QString intensityAverageFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/2. intensity pca on 130 meshes/mesh intensity pca average.txt");
	QString intensityEVectorFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/2. intensity pca on 130 meshes/mesh intensity pca eigenvectors.txt");
	
	QString combinedEVectorFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/3. combined pca on 130 meshes/mesh combined pca eigenvectors.txt");
	QString combinedEValuesFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/3. combined pca on 130 meshes/mesh combined pca eigenvalues.txt");
	QString wFileName = ("C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/3. combined pca on 130 meshes/w and nOfDataset.txt");
	
	QString combinedWeightsFileName = QFileDialog::getOpenFileName(this, "Load the weights", "C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/4. new instances", "*.txt");
	
	QString outputFolderName = QFileDialog::getExistingDirectory(this, "Load output folder", "C:/0.Data/2. Statistical model/2. mesh/b. 130 meshes/4. new instances");
	

	PCAMesh* PCA = new PCAMesh;
	PCA->SetNumberOfInstances(_nOfMeshInstances);
	PCA->SetNumberOfModes(_nOfMeshModes); 
	
	PCA->SetShapeAverageFileName(shapeAverageFileName );
	PCA->SetShapeEVectorsFileName(shapeEVectorFileName);
			
	PCA->SetIntensityAverageFileName(intensityAverageFileName);
	PCA->SetIntensityEVectorsFileName(intensityEVectorFileName);

	PCA->SetCombinedEVectorsFileName(combinedEVectorFileName);
	PCA->SetCombinedEValuesFileName(combinedEValuesFileName);
	PCA->SetCombinedWeightsFileName(combinedWeightsFileName);

	PCA->SetWFileName(wFileName);
	PCA->SetOutputFolder(outputFolderName);
	PCA->InstanceCreation();


}

/************************** FEM TOOL **************************/
// material properties
void MainWindow::loadMesh(){

	_meshFileNames = QFileDialog::getOpenFileNames(this, "Load meshes"," ","Abaqus quadratic volume meshes (*.inp);;Ansys quadratic volume meshes (*.cdb)");

}
void MainWindow::loadImage(){

	_imageFileNames = QFileDialog::getOpenFileNames(this, "Load images"," ","Images(*.mhd)");
}
void MainWindow::assignmentLawFirst(double value){

	_assignmentLawOne = value;
}
void MainWindow::assignmentLawSecond(double value){

	_assignmentLawTwo = value;
}
void MainWindow::assignToElements(){

	_nodeElementFlag = 1;
}
void MainWindow::assignToNodes(){

	_nodeElementFlag = 2;
}
void MainWindow::abaqusMatProp(){

	_fileType = 1;
}
void MainWindow::txtMatProp(){

	_fileType = 2;
}
void MainWindow::assign(){

	// define kind of assignment
	if (_nodeElementFlag != 1 && _nodeElementFlag != 2)
		std::cout << "define the kind of assignment" << std::endl;
	if (_fileType != 1 && _fileType != 2)
		std::cout << "define the file format" << std::endl;
		
	std::cout << "! be sure you are using calibrated or pseudo-calibrated images" << std::endl;
	std::cout << "all interpolated final grey values will be divided by 1000 - from [g/mm3 to mg/mm3]" << std::endl;


	// for each bone 
	double start = clock();

	for (int i=0; i<_meshFileNames.size(); i++){
		
		//std::cout << _meshFileNames[i].ascii() << std::endl; 
		//std::cout << _imageFileNames[i].ascii() << std::endl; 


		vtkPolyData* mesh = vtkPolyData::New();
		vtkDoubleArray* youngsModulus = vtkDoubleArray::New();

		// abaqus mesh
		if (_meshFileNames[i].endsWith(".inp")){
			MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
			meshReader->SetFileName(_meshFileNames[i]);
			std::cout << "mesh: " << _meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			mesh->DeepCopy(meshReader->GetOutput());
			delete meshReader;
		}
		
		// ansys mesh
		if (_meshFileNames[i].endsWith(".cdb")){ // Ansys quadratic volume meshes
			MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
			meshReader->SetFileName(_meshFileNames[i]);
			std::cout << "mesh: " << _meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			mesh->DeepCopy(meshReader->GetOutput());
			delete meshReader;		
		}
	
		/*
		// load the corresponding image (string check for stryker images)
		ImageHandler* imageHandler = new ImageHandler;
		QString temp = _meshFileNames[i];
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
				imageHandler->SetImageFileName(_imageFileNames[y]);
				std::cout << _imageFileNames[y].ascii() << std::endl;
			}
		}
		*/
		
		 
		ImageHandlerFloat* imageHandlerFloat = new ImageHandlerFloat;
	
		// load the correspondent image (string check for created instances)
		QString temp = _meshFileNames[i];
		if (temp.lastIndexOf("/") == -1){
				temp.remove("inp");
				temp.remove(0, temp.lastIndexOf("\\")+1);
				temp.replace(QString("."), QString("_"));
		}
		else {
			temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.replace(QString("."), QString("_"));
		}
		std::cout << "first " << temp.ascii() << std::endl;
		int length = temp.length();
		temp.remove(5,length);
		std::cout << "second " << temp.ascii() << std::endl;
		/*
		QString temp = _meshFileNames[i];
		if (temp.lastIndexOf("/") == -1){
			temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("\\")+1);
			temp.replace(QString("."), QString(""));
		}
		else {
			temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.replace(QString("."), QString(""));
		}
		*/

		
		for (int y=0; y<_imageFileNames.size(); y++){
			if (_imageFileNames[y].contains(temp)){
				imageHandlerFloat->SetImageFileName(_imageFileNames[y]);
				std::cout << "image: " << _imageFileNames[y].ascii() << std::endl;
			}
		}


		imageHandlerFloat->MetafileReaderUpdate();
		std::cout << "ITK to VTK to ITK coordinate system" << std::endl;
		imageHandlerFloat->ITKtoVTKtoITK(); // to have the same coordinate system
		
		std::cout << "extruding image" << std::endl;
		imageHandlerFloat->Extrusion(); // image extrusion
		
		
		// assign to the elements 

		if (_nodeElementFlag == 1){ //element

			FemAssignerElements* femAssigner = new FemAssignerElements; 
			// set mesh
			femAssigner->SetMesh(mesh);
			// set image
			femAssigner->SetImage(imageHandlerFloat->GetImage());
			// set laws
			femAssigner->SetAssignmentLawOne(_assignmentLawOne);
			femAssigner->SetAssignmentLawTwo(_assignmentLawTwo);
									
			// grey levels assignment
			std::cout << "assigning properties" << std::endl;
			femAssigner->GreyLevelAssignmentUpdate();
			// mechanical property assignment
			femAssigner->Update();
			youngsModulus->DeepCopy(femAssigner->GetYoungModulus());

			// get the apparent density
			vtkDoubleArray* appDensity = vtkDoubleArray::New();
			appDensity->DeepCopy(femAssigner->GetRhoApp());
			// save the apparent density
			QString temp = _meshFileNames[i];
			temp.replace(QString(".inp"), QString("_appDensity_elements.txt"));
			std::cout << temp.ascii() << std::endl;
			QFile outFile(temp);
			outFile.open(QIODevice::WriteOnly | QIODevice::Text);
			QTextStream writeFile(&outFile);

			for (int i=0; i<appDensity->GetNumberOfTuples();i++){
				if (i == appDensity->GetNumberOfTuples()-1)
					writeFile << appDensity->GetValue(i);
				else
					writeFile << appDensity->GetValue(i) << endl;
				
			}
			outFile.close();

			
			// extracting the first node of each element for rendering
			vtkPoints* firstNode = vtkPoints::New();
			vtkGenericCell* cell = vtkGenericCell::New();
			vtkIdList* idList;
			for (int i=0; i<mesh->GetNumberOfCells(); i++){
				// get the element
				mesh->GetCell(i, cell);
				idList = cell->GetPointIds();
				// get the first node
				double pt[3];
				mesh->GetPoint(idList->GetId(0), pt);
				firstNode->InsertNextPoint(pt);
			}

			
			// rendering with colorbar
			RenderingPointWithColorbar* rendering = new RenderingPointWithColorbar;
			rendering->SetPoints(firstNode);
			rendering->SetColorBar(femAssigner->GetYoungModulus());
			rendering->SetRenderer(_ren);
		//rendering->UpdateYoungModulus();

			delete femAssigner;
			delete rendering;
			std::cout << "properties assigned" << std::endl; 
		}

		
		// assigned to the nodes 
		else if (_nodeElementFlag == 2){ //nodes

			FemAssignerNodes* femAssigner = new FemAssignerNodes; 
			// set mesh	
			femAssigner->SetMesh(mesh);
			// set image
			femAssigner->SetImage(imageHandlerFloat->GetImage());
			// set laws
			femAssigner->SetAssignmentLawOne(_assignmentLawOne);
			femAssigner->SetAssignmentLawTwo(_assignmentLawTwo);
						
			// grey levels assignment
			std::cout << "assigning properties" << std::endl;
			femAssigner->GreyLevelAssignmentUpdate();
			// mechanical property assignment
			femAssigner->Update();
			youngsModulus->DeepCopy(femAssigner->GetYoungModulus());
			
			/*
			// get the ash density
			vtkDoubleArray* ashDensity = vtkDoubleArray::New();
			ashDensity->DeepCopy(femAssigner->GetRhoAsh());
			// save the apparent density
			QString temp = _meshFileNames[i];
			temp.replace(QString(".inp"), QString("_ashDensity_elements.txt"));
			std::cout << temp.ascii() << std::endl;
			QFile outFile(temp);
			outFile.open(QIODevice::WriteOnly | QIODevice::Text);
			QTextStream writeFile(&outFile);

			for (int i=0; i<ashDensity->GetNumberOfTuples();i++){
				if (i == ashDensity->GetNumberOfTuples()-1)
					writeFile << ashDensity->GetValue(i);
				else
					writeFile << ashDensity->GetValue(i) << endl;
				
			}
			outFile.close();
			*/
			/*
			// get the apparent density
			vtkDoubleArray* appDensity = vtkDoubleArray::New();
			appDensity->DeepCopy(femAssigner->GetRhoApp());
			// save the apparent density
			QString temp = _meshFileNames[i];
			temp.replace(QString(".inp"), QString("_appDensity_elements.txt"));
			std::cout << temp.ascii() << std::endl;
			QFile outFile(temp);
			outFile.open(QIODevice::WriteOnly | QIODevice::Text);
			QTextStream writeFile(&outFile);

			for (int i=0; i<appDensity->GetNumberOfTuples();i++){
				if (i == appDensity->GetNumberOfTuples()-1)
					writeFile << appDensity->GetValue(i);
				else
					writeFile << appDensity->GetValue(i) << endl;
				
			}
			outFile.close();
			*/


			// rendering with colorbar
			RenderingPointWithColorbar* rendering = new RenderingPointWithColorbar;
			rendering->SetPoints(mesh->GetPoints());
			rendering->SetColorBar(femAssigner->GetYoungModulus());
			rendering->SetRenderer(_ren);
		//rendering->UpdateYoungModulus();

			delete femAssigner;
			delete rendering;
			std::cout << "properties assigned" << std::endl; 
		}
		
		if (_fileType == 1){ // save as abaqus

			std::cout << "option disabled" << std::endl;

			/*
			MeshWriterAbaqus* writer = new MeshWriterAbaqus;
			if (_nodeElementFlag == 1){ //element
				if (_meshFileNames[i].contains(".inp"))
					_meshFileNames[i].replace(QString(".inp"), QString("_mechProp_elements.inp"));
				else if (_meshFileNames[i].contains(".cdb"))
					_meshFileNames[i].replace(QString(".cdb"), QString("_mechProp_elements.inp"));
			}
			else if (_nodeElementFlag == 2){ //nodes
				if (_meshFileNames[i].contains(".inp"))
					_meshFileNames[i].replace(QString(".inp"), QString("_mechProp_nodes.inp"));
				else if (_meshFileNames[i].contains(".cdb"))
					_meshFileNames[i].replace(QString(".cdb"), QString("_mechProp_nodes.inp"));
			}
			writer->SetFileName(_meshFileNames[i]);
			writer->MeshOn();
			writer->SetMesh(mesh);
			if (_nodeElementFlag == 1)
				writer->YoungsModulusElemOn();
			else if (_nodeElementFlag == 2)
				writer->YoungsModulusNodeOn();
			writer->SetYoungsModulus(youngsModulus);
			writer->Update();

			delete writer;
			*/
		
		}
		else if (_fileType == 2){ // save as txt
			
			if (_nodeElementFlag == 1)
				_meshFileNames[i].replace(QString(".inp"), QString("_mechProp_elements.txt"));
			else if (_nodeElementFlag == 2)
				_meshFileNames[i].replace(QString(".inp"), QString("_mechProp_nodes.txt"));
						
			//temp.append("_mechProp_nodes.txt"); // original/generalization/representation instances
			//std:cout << "output file: " << temp.toAscii().data() << std::endl;
			//QFile outFile(temp);

			QFile outFile(_meshFileNames[i]);
			outFile.open(QIODevice::WriteOnly | QIODevice::Text);
			QTextStream writeFile(&outFile);

			for (int i=0; i<youngsModulus->GetNumberOfTuples();i++){
				if (i == youngsModulus->GetNumberOfTuples()-1)
					writeFile << youngsModulus->GetValue(i);
				else
					writeFile << youngsModulus->GetValue(i) << endl;
				
			}
			outFile.close();

		}
		
		// cleaning
		delete imageHandlerFloat;

		std::cout << std::endl;
		std::cout << "file written" << std::endl;
		
	
	}
	double end = clock();
	double total = (end-start)/CLOCKS_PER_SEC;
	std::cout << "Computation time: " << total << " sec" << std::endl;
}

// boundary conditions
void MainWindow::pickApplicationPoints(){

	_fileName = QFileDialog::getOpenFileName(this, "Load bone meshes","C:/0.Data/2011_FemaleMale_Study/6. new instances/female/0.reference","Abaqus quadratic volume meshes (*.inp);;Ansys quadratic volume meshes(*.cdb)");
	
	vtkPolyData* boneSurface = vtkPolyData::New();
		
		// abaqus
		if (_fileName.endsWith(".inp")){
			MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
			meshReader->SetFileName(_fileName);
			std::cout << _fileName.toAscii().data() << std::endl; 
			meshReader->Update();
			boneSurface->DeepCopy(meshReader->GetOutput());		
		}
		// ansys
		else if (_fileName.endsWith(".cdb")){
			MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
			meshReader->SetFileName(_fileName);
			std::cout << _fileName.toAscii().data() << std::endl; 
			meshReader->Update();
			boneSurface->DeepCopy(meshReader->GetOutput());	
		}

		// extracting outer surface
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		surfaceExtractor->SetVolumeMesh(boneSurface);
		surfaceExtractor->Update();
		
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(surfaceExtractor->GetSurfaceMesh());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();
	
	
	// picker
	std::cout << "position the mouse on the point + press p" << std::endl;
	_pointPicker->SetMesh(boneSurface);
	_pointPicker->SetRenderer(_ren);
	double color[3];
	color[0] = 1.0; color[1] = 0.0; color[2] = 0.0;
	_pointPicker->SetColor(color);
	_pointPicker->SetPicker(_cellPicker);
	_coords = vtkPoints::New();
	_pointPicker->SetPoints(_coords);
	_pointPicker->Update();

}
void MainWindow::saveApplicationPoints(){

	_cellPicker->RemoveAllObservers();
	
	// file name 
	if (_fileName.endsWith(".inp"))
		_fileName.replace(QString(".inp"), QString("_points.txt"));
	else if (_fileName.endsWith(".cdb"))
		_fileName.replace(QString(".cdb"), QString("_points.txt"));
	
	// print out
	for (int i=0; i<_coords->GetNumberOfPoints(); i++){
		double coord[3];
		_coords->GetPoint(i, coord);
		std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
	}
	
	// save points
	PointWriterXyz* pointWriter = new PointWriterXyz;
	pointWriter->SetInput(_coords);
	pointWriter->SetFileName(_fileName.toAscii().data());
	pointWriter->Update();
	std::cout << ".txt load file written" << std::endl;

}
void MainWindow::propagateApplicationPoints(){

	/*
	std::cout << "! string check done on instance_n.inp and svf_n.mhd" << std::endl;
	
	QString pointFileName = QFileDialog::getOpenFileName(this, "Load points","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/0.reference","*.txt");
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load meshes","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/4.meshes","*.inp");
	QStringList SVFfileNames = QFileDialog::getOpenFileNames(this, "Load SVF","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/1.inverted svf","*.mhd");
	//QString pointFileName = QFileDialog::getOpenFileName(this, "Load points","C:/0.Data/test data/fe simulation/reference","*.txt");
	//QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load meshes","C:/0.Data/test data/fe simulation","*.inp");
	//QStringList SVFfileNames = QFileDialog::getOpenFileNames(this, "Load SVF","C:/0.Data/test data/fe simulation","*.mhd");
	

	// for each bone 
	for (int i=0; i<meshFileNames.size(); i++){
		
		// load points to propagate
		PointReaderXyz* pointReader = new PointReaderXyz;
		pointReader->SetFileName(pointFileName);
		std::cout << "application points: " << pointFileName.toAscii().data() << std::endl;
		pointReader->Update();

		FemForce* propagatePoints = new FemForce;
		propagatePoints->SetPoints(pointReader->GetOutput());
		for (int a=0; a<propagatePoints->GetPoints()->GetNumberOfPoints(); a++){
			double pt[3]; 
			propagatePoints->GetPoints()->GetPoint(a,pt);
			std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
		}

		// load mesh
		MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
		meshReader->SetFileName(meshFileNames[i]);
		std::cout << "mesh: " << meshFileNames[i].ascii() << std::endl; 
		meshReader->Update();
		propagatePoints->SetMesh(meshReader->GetOutput());

		// load the correspondent SVF
		ImageHandler* fieldHandler = new ImageHandler;
		QString temp = meshFileNames[i];
		int length = temp.length();
		if (temp.lastIndexOf("/") == -1){
			temp.remove(".inp");
			temp.remove(0, temp.lastIndexOf("\\")+1);
			temp.remove(5, length);
		}
		else {
			temp.remove("inp");
			temp.remove(".inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.remove(5, length);
		}

		for (int y=0; y<SVFfileNames.size(); y++){
			if (SVFfileNames[y].contains(temp)){
				fieldHandler->SetFieldFileName(SVFfileNames[y]);
				std::cout << "svf: " << SVFfileNames[i].ascii() << std::endl;
			}
		}
		fieldHandler->FieldReaderUpdate();
		propagatePoints->SetSVF(fieldHandler->GetField());
		
		// point propagations
		propagatePoints->PropagateApplicationPoints();

		// save propagated points
		if (meshFileNames[i].endsWith(".inp"))
			meshFileNames[i].replace(QString(".inp"), QString("_applicationPoints_neck.txt"));
		else if (meshFileNames[i].endsWith(".cdb"))
			meshFileNames[i].replace(QString(".cdb"), QString("_applicationPoints_neck.txt"));
		else 
			meshFileNames[i].append("test.txt");
		PointWriterXyz* pointWriter = new PointWriterXyz;
		pointWriter->SetFileName(meshFileNames[i]);
		pointWriter->SetInput(propagatePoints->GetPoints());
		pointWriter->Update();
		
		
		// save propagated points with index
		//meshFileNames[i].replace(QString("_propagatedPt.txt"), QString("_propagatedPt_index.txt"));
		if (meshFileNames[i].endsWith(".inp"))
			meshFileNames[i].replace(QString(".inp"), QString("_applicationPoints_neck.txt"));
		else if (meshFileNames[i].endsWith(".cdb"))
			meshFileNames[i].replace(QString(".cdb"), QString("_applicationPoints_neck.txt"));
		PointWriterXyzId* pointIDWriter = new PointWriterXyzId;
		pointIDWriter->SetFileName(meshFileNames[i]);
		pointIDWriter->SetInput(propagatePoints->GetPoints());
		pointIDWriter->SetIdVector(propagatePoints->GetPointsID());
		pointIDWriter->Update();

		std::cout << std::endl; 

		// cleaning
		
		delete propagatePoints;
		delete meshReader;
		delete fieldHandler;
		//delete pointWriter;
		delete pointReader;
		delete pointIDWriter;
	}
	std::cout << "propagation done" << std::endl; 
	*/


	/*
	///////////////////////////////////////////////////////////////////
	///////////////////// get mesh closest points /////////////////////
	///////////////////////////////////////////////////////////////////

	QString pointFileName = QFileDialog::getOpenFileName(this, "Load points"," ","*.txt");
	QString meshFileName = QFileDialog::getOpenFileName(this, "Load meshes"," ","*.inp");
	//QString meshFileName = QFileDialog::getOpenFileName(this, "Load meshes"," ","*.cdb");
	
	// load mesh
	MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
	//MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
	meshReader->SetFileName(meshFileName);
	std::cout << "mesh: " << meshFileName.ascii() << std::endl; 
	meshReader->Update();
	MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
	surfaceExtractor->SetVolumeMesh(meshReader->GetOutput());
	surfaceExtractor->Update();
	
	// load points to propagate
	PointReaderXyz* pointReader = new PointReaderXyz;
	pointReader->SetFileName(pointFileName);
	std::cout << "application points: " << pointFileName.toAscii().data() << std::endl;
	pointReader->Update();
	
	// propagate points
	FemForce* propagatePoints = new FemForce;
	propagatePoints->SetMesh(surfaceExtractor->GetSurfaceMesh());
	propagatePoints->SetPoints(pointReader->GetOutput());

	vtkPoints* points = vtkPoints::New();
		 
	for (int a=0; a<propagatePoints->GetPoints()->GetNumberOfPoints(); a++){
		
		double pt[3]; 
		propagatePoints->GetPoints()->GetPoint(a,pt);
		std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;

		int pointIndex;
		double closestX; double closestY; double closestZ;
		propagatePoints->findClosestNode (surfaceExtractor->GetSurfaceMesh(), pt, pointIndex, closestX, closestY, closestZ);
		std::cout << "closest point" << std::endl;
		double coord[3];
		coord[0] = closestX;
		coord[1] = closestY;
		coord[2] = closestZ;
		
		std::cout << closestX << ' ' << closestY << ' ' << closestZ << std::endl;
		std::cout << pointIndex << std::endl;
			
		// put it back in vtkPoints
		points->InsertPoint(a,coord);

	}
	*/
	
	
	// extract points with same ID from all the meshes //
	
	QString pointFileName = QFileDialog::getOpenFileName(this, "Load points"," ","*.txt");
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load meshes"," ","*.inp");
	//QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load meshes"," ","*.cdb");

	// load points and get IDs
	PointReaderXyzId* pointReader = new PointReaderXyzId;
	pointReader->SetFileName(pointFileName);
	std::cout << "application points: " << pointFileName.toAscii().data() << std::endl;
	pointReader->Update();
	vtkDoubleArray* IDs = vtkDoubleArray::New();
	IDs->DeepCopy(pointReader->GetIdVector());
	
	for (int i=0; i<meshFileNames.size(); i++){

		// load mesh and extract nodes
		MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
		//MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
		meshReader->SetFileName(meshFileNames[i]);
		std::cout << "mesh: " << meshFileNames[i].ascii() << std::endl; 
		meshReader->Update();
		vtkPoints* nodes = vtkPoints::New();
		nodes = meshReader->GetOutput()->GetPoints();
	
		// get nodes
		vtkPoints* points = vtkPoints::New();
		for (int a=0; a<IDs->GetNumberOfTuples(); a++){
			int id = IDs->GetValue(a);
			id=id-1;
			double pt[3];
			nodes->GetPoint(id, pt);
			points->InsertNextPoint(pt);
		}

		// write points 
		QString outputFileName = meshFileNames[i];
		outputFileName.replace(".inp", "_applicationPoints_neck.txt");
		std::cout << outputFileName.toAscii().data() << std::endl;
		PointWriterXyzId* pointWriter = new PointWriterXyzId;
		pointWriter->SetInput(points);
		pointWriter->SetIdVector(IDs);
		pointWriter->SetFileName(outputFileName);
		pointWriter->Update();
	}
	
	

}
void MainWindow::walking(){

	QStringList pointsFileNames = QFileDialog::getOpenFileNames(this, "Load the propagated application points for walking","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/6.applPoints","*.txt");
	
	// for each bone 
	for (int i=0; i<pointsFileNames.size(); i++){
		
		// load propagated points
		PointReaderXyzId* pointReader = new PointReaderXyzId;
		pointReader->SetFileName(pointsFileNames[i]);
		std::cout << "points: " << pointsFileNames[i].toAscii().data() << std::endl;
		pointReader->Update();
		
		// walking forces (hip contanct + muscles)
		FemForce* walk = new FemForce;
		walk->SetPoints(pointReader->GetOutput());
		walk->SetPointsId(pointReader->GetIdVector());
		walk->WalkingFromFemulLenght(); // walking
		//walk->Standing(); // standing (comparison paper = fem validation)

		// writing forces 
		PointWriterXyzId* pointWriter = new PointWriterXyzId;
		pointWriter->SetInput(walk->GetPoints());
		pointWriter->SetIdVector(walk->GetPointsID());
		QString temp = pointsFileNames[i];
		if (temp.lastIndexOf("/") == -1){
				//temp.remove("inp");
				temp.remove(0, temp.lastIndexOf("\\")+1);
				temp.replace(QString("."), QString("_"));
		}
		else {
			//temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.replace(QString("."), QString("_"));
		}
		int length = temp.length();
		temp.remove(5,length);
		temp.append("_forcesBC.txt");
		//pointsFileNames[i].replace(QString("applicationPoints"), QString("forcesBC"));
		std::cout << "output: " << temp.toAscii().data() << std::endl;
		pointWriter->SetFileName(temp);
		pointWriter->Update();

		// cleaning
		delete pointReader;
		delete walk;
		delete pointWriter;
	}

	std::cout << "walking forces and boundary conditions calculated" << std::endl;

}
void MainWindow::falling(){


	QStringList pointsFileNames = QFileDialog::getOpenFileNames(this, "Load the propagated application points for falling","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/6.applPoints","*.txt");
	QStringList loadFileNames = QFileDialog::getOpenFileNames(this, "Load the forces for standing (bodyweight)","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/7.forces","*.txt");
	
	// for each bone 
	for (int i=0; i<pointsFileNames.size(); i++){
		
		// load propagated points
		PointReaderXyzId* pointReader = new PointReaderXyzId;
		pointReader->SetFileName(pointsFileNames[i]);
		std::cout << "points: " << pointsFileNames[i].toAscii().data() << std::endl;
		pointReader->Update();

		// load body weight
		PointReaderXyzId* weightPointReader = new PointReaderXyzId;
		QString temp = pointsFileNames[i];
		temp.replace("","");
		if (temp.lastIndexOf("/") == -1){
				temp.remove(0, temp.lastIndexOf("\\")+1);
				temp.replace(QString("applicationPoints"), QString("forcesBC"));
		}
		else {
			temp.remove(0, temp.lastIndexOf("\\")+1);
			temp.replace(QString("applicationPoints"), QString("forcesBC"));
		}
		
		for (int y=0; y<loadFileNames.size(); y++){
			if (loadFileNames[y].contains(temp)){
				std::cout << loadFileNames[i].ascii() << std::endl;
				weightPointReader->SetFileName(loadFileNames[y]);
			}
		}
		weightPointReader->Update();
		
		// falling forces
		FemForce* fall = new FemForce;
		fall->SetPoints(pointReader->GetOutput());
		fall->SetPointsId(pointReader->GetIdVector());
		
		double bodyWeight = weightPointReader->GetIdVector()->GetValue(9);
		std::cout << bodyWeight << std::endl;
		fall->SetBodyWeight(bodyWeight);
		
		fall->SideFalling();

		// writing forces 
		PointWriterXyzId* pointWriter = new PointWriterXyzId;
		pointWriter->SetInput(fall->GetPoints());
		pointWriter->SetIdVector(fall->GetPointsID());
		pointsFileNames[i].replace(QString("applicationPoints"), QString("forcesBC_falling"));
		std::cout << pointsFileNames[i].toAscii().data() << std::endl;
		pointWriter->SetFileName(pointsFileNames[i]);
		pointWriter->Update();

		// cleaning
		delete pointReader;
		delete weightPointReader;
		delete fall;
		delete pointWriter;
	}

	std::cout << "falling forces and boundary conditions calculated" << std::endl;

}

// finite element simulation
void MainWindow::createAbaqusInp(){

	
	// loading
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load bone mesh", "C:/0.Data/2011_FemaleMale_Study/6. new instances/male/4.meshes","Abaqus quadratic volume meshes (*.inp)");
	QStringList propertiesFileNames = QFileDialog::getOpenFileNames(this, "Load mechanical properties file","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/5.nodeMechProp","*.txt");
	QStringList loadFileNames = QFileDialog::getOpenFileNames(this, "Load force and boundary condition file","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/7.forces","*.txt");
	QStringList appFileNames = QFileDialog::getOpenFileNames(this, "Load application point file with ID","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/6.applPoints","*.txt");
	
	// ==============================================================================
	// for comparison paper - fem validation (the stance code is below commented out)
	
	for (int i=0; i<meshFileNames.size(); i++){
	
		FemAbaqusInpWriter* inpWriter = new FemAbaqusInpWriter;
	
		// mesh
		std::cout << "mesh" << std::endl; 
		MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
		//std::cout << meshFileNames.ascii() << std::endl; 
		std::cout << meshFileNames[i].toAscii().data() << std::endl; 
		meshReader->SetFileName(meshFileNames[i]);
		//meshReader->SetFileName(meshFileNames);
		meshReader->Update();
		inpWriter->SetMesh(meshReader->GetOutput());
	
		// string in fileNames
		QString temp = meshFileNames[i];
		if (temp.lastIndexOf("/") == -1){
				temp.remove("inp");
				temp.remove(0, temp.lastIndexOf("\\")+1);
				temp.replace(QString("."), QString("_"));
		}
		else {
			temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.replace(QString("."), QString("_"));
		}
		int length = temp.length();
		temp.remove(5,length);
		std::cout << "------ temp: " << temp.toAscii().data() << std::endl;

		// mechanical properties
		std::cout << "mechanical properties" << std::endl; 
		vtkDoubleArray* mechProp = vtkDoubleArray::New(); 
		for (int y=0; y<propertiesFileNames.size(); y++){
			if (propertiesFileNames[y].contains(temp)){
				std::cout << propertiesFileNames[y].toAscii().data() << std::endl; 
				QFile inFile(propertiesFileNames[y]);
				//std::cout << propertiesFileNames.ascii() << std::endl; 
				//QFile inFile(propertiesFileNames);
				inFile.open(QIODevice::ReadOnly | QIODevice::Text);
				QTextStream file(&inFile);
				while (!file.atEnd()){
					double temp;
					file >> temp;
					mechProp->InsertNextValue(temp);
				}
				inFile.close();
				inpWriter->SetMechProp(mechProp);
			}
		}
		
		
		// boundary conditions and load
		std::cout << "force and boundary conditions" << std::endl; 
		PointReaderXyzId* pointReader = new PointReaderXyzId;
		for (int y=0; y<loadFileNames.size(); y++){
			if (loadFileNames[y].contains(temp)){
				std::cout << loadFileNames[y].toAscii().data() << std::endl;
				pointReader->SetFileName(loadFileNames[y]);
			}
		}
		//std::cout << loadFileNames.ascii() << std::endl;
		//pointReader->SetFileName(loadFileNames);
		pointReader->Update();
		vtkPoints* coords = vtkPoints::New();
		coords = pointReader->GetOutput();
		vtkDoubleArray* id = vtkDoubleArray::New(); 
		id = pointReader->GetIdVector(); 
		inpWriter->SetBoundaryConditions(coords);
		inpWriter->SetBoundaryConditionsId(id);

		//application points
		std::cout << "force and boundary condition application points" << std::endl;
		PointReaderXyzId* appPointReader = new PointReaderXyzId;
		for (int y=0; y<appFileNames.size(); y++){
			if (appFileNames[y].contains(temp)){
				std::cout << appFileNames[y].toAscii().data() << std::endl;
				appPointReader->SetFileName(appFileNames[y]);
			}
		}
		//std::cout << appFileNames.ascii() << std::endl;
		//appPointReader->SetFileName(appFileNames);
		appPointReader->Update();
		inpWriter->SetApplicationPointsID(appPointReader->GetIdVector());
		
		// write file
		std::cout << "writing abaqus input file" << std::endl;	
		//meshFileNames[i].replace(QString(".inp"), QString("_abaqus.inp"));
		//meshFileNames.replace(QString(".inp"), QString("_abaqus.inp"));
		//inpWriter->SetFileName(meshFileNames[i]);
		//inpWriter->SetFileName(meshFileNames);
		temp.append("_abaqus.inp");
		inpWriter->SetFileName(temp);
		//inpWriter->WriteNodeWalkingInp(); 
		inpWriter->WriteNodeStandingInp(); // for comparison paper - fem validation

		std::cout << "output file name: " << temp.ascii() << std::endl;
		std::cout << "abaqus input file written" << std::endl;

		// cleaning
		delete meshReader;
		delete pointReader;
		delete appPointReader;
		delete inpWriter;
	}
	
	//======================================================================
/*

// for stance
	for (int i=0; i<meshFileNames.size(); i++){
	
		FemAbaqusInpWriter* inpWriter = new FemAbaqusInpWriter;
	
		// mesh
		std::cout << "mesh" << std::endl; 
		MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
		std::cout << meshFileNames[i].toAscii().data() << std::endl; 
		meshReader->SetFileName(meshFileNames[i]);
		meshReader->Update();
		inpWriter->SetMesh(meshReader->GetOutput());
	
		// string in fileNames
		QString temp = meshFileNames[i];
		if (temp.lastIndexOf("/") == -1){
				temp.remove("inp");
				temp.remove(0, temp.lastIndexOf("\\")+1);
				temp.replace(QString("."), QString("_"));
		}
		else {
			temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.replace(QString("."), QString("_"));
		}
		
		// mechanical properties
		std::cout << "mechanical properties" << std::endl; 
		vtkDoubleArray* mechProp = vtkDoubleArray::New(); 
		for (int y=0; y<propertiesFileNames.size(); y++){
			if (propertiesFileNames[y].contains(temp)){
				std::cout << propertiesFileNames[y].toAscii().data() << std::endl; 
				QFile inFile(propertiesFileNames[y]);
				inFile.open(QIODevice::ReadOnly | QIODevice::Text);
				QTextStream file(&inFile);
				while (!file.atEnd()){
					double temp;
					file >> temp;
					mechProp->InsertNextValue(temp);
				}
				inFile.close();
				inpWriter->SetMechProp(mechProp);
			}
		}
		
		
		// boundary conditions and load
		std::cout << "force and boundary conditions" << std::endl; 
		PointReaderXyzId* pointReader = new PointReaderXyzId;
		for (int y=0; y<loadFileNames.size(); y++){
			if (loadFileNames[y].contains(temp)){
				std::cout << loadFileNames[y].toAscii().data() << std::endl;
				pointReader->SetFileName(loadFileNames[y]);
			}
		}
		pointReader->Update();
		vtkPoints* coords = vtkPoints::New();
		coords = pointReader->GetOutput();
		vtkDoubleArray* id = vtkDoubleArray::New(); 
		id = pointReader->GetIdVector(); 
		inpWriter->SetBoundaryConditions(coords);
		inpWriter->SetBoundaryConditionsId(id);

		//application points
		std::cout << "force and boundary condition application points" << std::endl;
		PointReaderXyzId* appPointReader = new PointReaderXyzId;
		for (int y=0; y<appFileNames.size(); y++){
			if (appFileNames[y].contains(temp)){
				std::cout << appFileNames[y].toAscii().data() << std::endl;
				appPointReader->SetFileName(appFileNames[y]);
			}
		}
		appPointReader->Update();
		inpWriter->SetApplicationPointsID(appPointReader->GetIdVector());
		
		// write file
		std::cout << "writing abaqus input file" << std::endl;	
		meshFileNames[i].replace(QString(".inp"), QString("_abaqus_m_500_stance.inp"));
		inpWriter->SetFileName(meshFileNames[i]);
		inpWriter->WriteNodeWalkingInp();
		std::cout << "abaqus input file written" << std::endl;

		// cleaning
		delete meshReader;
		delete pointReader;
		delete appPointReader;
		delete inpWriter;
	}
	*/
	
	
}
void MainWindow::createAbaqusInpFalling(){

	// loading 
	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load bone mesh", "C:/0.Data/2011_FemaleMale_Study/6. new instances/male/4.meshes","Abaqus quadratic volume meshes (*.inp)");
	QStringList propertiesFileNames = QFileDialog::getOpenFileNames(this, "Load mechanical properties file","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/5.nodeMechProp","*.txt");
	QStringList loadFileNames = QFileDialog::getOpenFileNames(this, "Load force and boundary condition file","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/7.forces/falling","*.txt");
	QStringList appFileNames = QFileDialog::getOpenFileNames(this, "Load application point file with ID","C:/0.Data/2011_FemaleMale_Study/6. new instances/male/6.applPoints","*.txt");
	
	
	for (int i=0; i<meshFileNames.size(); i++){
	
		FemAbaqusInpWriter* inpWriter = new FemAbaqusInpWriter;
	
		// mesh
		std::cout << "mesh" << std::endl; 
		MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
		std::cout << meshFileNames[i].toAscii().data() << std::endl; 
		meshReader->SetFileName(meshFileNames[i]);
		meshReader->Update();
		inpWriter->SetMesh(meshReader->GetOutput());
	
		// string in fileNames
		QString temp = meshFileNames[i];
		if (temp.lastIndexOf("/") == -1){
				temp.remove("inp");
				temp.remove(0, temp.lastIndexOf("\\")+1);
				temp.replace(QString("."), QString("_"));
		}
		else {
			temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.replace(QString("."), QString("_"));
		}
		
		// mechanical properties
		std::cout << "mechanical properties" << std::endl; 
		vtkDoubleArray* mechProp = vtkDoubleArray::New(); 
		for (int y=0; y<propertiesFileNames.size(); y++){
			if (propertiesFileNames[y].contains(temp)){
				std::cout << propertiesFileNames[y].toAscii().data() << std::endl; 
				QFile inFile(propertiesFileNames[y]);
				inFile.open(QIODevice::ReadOnly | QIODevice::Text);
				QTextStream file(&inFile);
				while (!file.atEnd()){
					double temp;
					file >> temp;
					mechProp->InsertNextValue(temp);
				}
				inFile.close();
				inpWriter->SetMechProp(mechProp);
			}
		}
		
		
		// boundary conditions and load
		std::cout << "force and boundary conditions" << std::endl; 
		PointReaderXyzId* pointReader = new PointReaderXyzId;
		for (int y=0; y<loadFileNames.size(); y++){
			if (loadFileNames[y].contains(temp)){
				std::cout << loadFileNames[y].toAscii().data() << std::endl;
				pointReader->SetFileName(loadFileNames[y]);
			}
		}
		pointReader->Update();
		vtkPoints* coords = vtkPoints::New();
		coords = pointReader->GetOutput();
		vtkDoubleArray* id = vtkDoubleArray::New(); 
		id = pointReader->GetIdVector(); 
		inpWriter->SetBoundaryConditions(coords);
		inpWriter->SetBoundaryConditionsId(id);

		//application points
		std::cout << "force and boundary condition application points" << std::endl;
		PointReaderXyzId* appPointReader = new PointReaderXyzId;
		for (int y=0; y<appFileNames.size(); y++){
			if (appFileNames[y].contains(temp)){
				std::cout << appFileNames[y].toAscii().data() << std::endl;
				appPointReader->SetFileName(appFileNames[y]);
			}
		}
		appPointReader->Update();
		inpWriter->SetApplicationPointsID(appPointReader->GetIdVector());
		
		// write file
		std::cout << "writing abaqus input file" << std::endl;	
		meshFileNames[i].replace(QString(".inp"), QString("_abaqus_m_500_falling.inp"));
		inpWriter->SetFileName(meshFileNames[i]);
		inpWriter->WriteNodeFallingInp();
		std::cout << "abaqus input file written" << std::endl;

		// cleaning
		delete meshReader;
		delete pointReader;
		delete appPointReader;
		delete inpWriter;

		std::cout << std::endl; 

	}
}
void MainWindow::runSimulation(){

	// in cmd: 
	// cd "folder for output"
	// abaqus j=jobName.inp interactive

	// abaqus output files will go to the release or debug forder
			
		
	_fileNames = QFileDialog::getOpenFileNames(this, "Load abaqus input files","C:/0.Data/test data/abaqus file creation","*.inp");
	
	for (int i=0; i<_fileNames.size(); i++){
		
		std::cout << _fileNames[i].toAscii().data() << std::endl;
		QString temp = _fileNames[i];
		if (temp.lastIndexOf("/") == -1){
				temp.remove("inp");
				temp.remove(0, temp.lastIndexOf("\\")+1);
				temp.replace(QString("."), QString("_"));
		}
		else {
			temp.remove("inp");
			temp.remove(0, temp.lastIndexOf("/")+1);
			temp.replace(QString("."), QString("_"));
		}

		
		int length = temp.length();
		//temp.replace(QString("_abaqus_"), QString("")); // standard
		temp.replace(QString("abaqus_"), QString("")); // for ben
		//temp.remove("_abaqus_");
		

		QString arguments;
		arguments = ("\"C:/SIMULIA/Abaqus/Commands/abq6101.bat\" j=job_");
		arguments.append(temp);
		arguments.append(" inp=");
		arguments.append("\"");
		arguments.append(_fileNames[i]);
		arguments.append("\"");
		arguments.replace(QString("\\"), QString("/"));
		arguments.append(" interactive");
		
		std::cout << "to abaqus :" << arguments.toAscii().data() << std::endl;
		
		double start = clock();
		QProcess myProcess;
		myProcess.start(arguments.toAscii().data());
		myProcess.waitForFinished(2000000);
		double end = clock();
		double total = (end-start)/CLOCKS_PER_SEC;
		std::cout << "Computation time: " << total << " sec. (" << total/60 << " min.)" << std::endl;

		std::cout << "simulation done" << std::endl;
		
		
	}
}
void MainWindow::extractNeckIntensities(){


	QString referenceMaskFileName = QFileDialog::getOpenFileName(this, "Load reference mask","C:/0.Data/2011_FemaleMale_Study/6. new instances/female/0.reference","Images(*.mhd)");
	
	QStringList svfFileNames = QFileDialog::getOpenFileNames(this, "Load inverted svf","C:/0.Data/2011_FemaleMale_Study/6. new instances/female/1.inverted svf","(*.mhd)");
	QStringList imageFileNames = QFileDialog::getOpenFileNames(this, "Load images","C:/0.Data/2011_FemaleMale_Study/6. new instances/female/3.pseudocalibrated","Images(*.mhd)");
	
	// reference mask
	ImageHandlerFloat* referenceMask = new ImageHandlerFloat;
	referenceMask->SetImageFileName(referenceMaskFileName);
	referenceMask->MetafileReaderUpdate();

	for (int i=0; i<svfFileNames.size(); i++){
		
		
		ImageHandlerFloat* imageHandlerFloat = new ImageHandlerFloat;
		
		// reference
		imageHandlerFloat->SetMask(referenceMask->GetImage());
		
		// svf
		VectorImageHandler* svf = new VectorImageHandler;
		std::cout << svfFileNames[i].ascii() << std::endl;
		svf->SetFieldFileName(svfFileNames[i]);
		svf->MetafileReaderUpdate();
		std::cout << "inverting the svf and tranforming to dvf" << std::endl;
		svf->VFtoDVFinverted();
		imageHandlerFloat->SetSVF(svf->GetField());

		// correspondent image
		QString temp = svfFileNames[i];
		temp.remove(".mhd");
		temp.remove(0, temp.lastIndexOf("_"));
		temp.append("_pseudoCalibrated");
		std::cout << temp.ascii() << std::endl;
		for (int y=0; y<imageFileNames.size(); y++){
			if (imageFileNames[y].contains(temp)){
				std::cout << imageFileNames[y].ascii() << std::endl;
				imageHandlerFloat->SetImageFileName(imageFileNames[y]);
			}
		}
		imageHandlerFloat->MetafileReaderUpdate();

		// extracting the neck
		std::cout << "extracting the neck volume" << std::endl; 
		imageHandlerFloat->ExtractVolumeImage();

		// saving the gray level vector
		vnl_vector<double> neckIntensities;
		neckIntensities = imageHandlerFloat->GetVector();

		VnlWriterVector* writer = new VnlWriterVector;
		svfFileNames[i].replace("inverted_svf", "instance");
		svfFileNames[i].replace(".mhd", "_neckIntensities.txt");
		std::cout << svfFileNames[i].ascii() << std::endl;
		writer->SetFileName(svfFileNames[i]);
		writer->SetVnlVector(neckIntensities);
		writer->Update();

		std::cout << std::endl;

		delete imageHandlerFloat;
		delete svf;
		delete writer;
	}

	delete referenceMask;

		

}

// implant fitting
void MainWindow::loadBoneMesh(){

	_meshFileNames = QFileDialog::getOpenFileNames(this, "Load bone meshes","C:/0.Data/test data/bone implant fitting","STL meshes(*.stl);;Abaqus quadratic volume meshes (*.inp);;Ansys quadratic volume meshes(*.cdb)");
}
void MainWindow::loadImplantMesh(){

	_fileName = QFileDialog::getOpenFileName(this, "Load implant mesh", "C:/0.Data/test data/bone implant fitting", "*.stl");

}
void MainWindow::proximal(){

	_positionFlag = 1;
}
void MainWindow::diaphyseal(){

	_positionFlag = 2;
}
void MainWindow::distal(){

	_positionFlag = 3;
}
void MainWindow::computeAndSave(){

	// read implant mesh
	vtkSTLReader* implantReader = vtkSTLReader::New();
	std::cout << _fileName.toAscii().data() << std::endl; 
	implantReader->SetFileName(_fileName.toAscii().data());
	implantReader->Update();
	std::cout << "nodes: " << implantReader->GetOutput()->GetNumberOfPoints() << std::endl;
	std::cout << "elements: " <<  implantReader->GetOutput()->GetNumberOfCells() << std::endl;
	
	// for each bone
	for (int i=0; i<_meshFileNames.size(); i++){
		
		// read bone mesh
		vtkPolyData* boneSurface = vtkPolyData::New();
		
		// abaqus
		if (_meshFileNames[0].endsWith(".inp")){
			MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
			meshReader->SetFileName(_meshFileNames[i]);
			std::cout << _meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			boneSurface->DeepCopy(meshReader->GetOutput());		
		}
		// ansys
		else if (_meshFileNames[0].endsWith(".cdb")){
			MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
			meshReader->SetFileName(_meshFileNames[i]);
			std::cout << _meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			boneSurface->DeepCopy(meshReader->GetOutput());	
		}
		// stl
		else if (_meshFileNames[0].endsWith(".stl")){
			vtkSTLReader* meshReader = vtkSTLReader::New();
			meshReader->SetFileName(_meshFileNames[i]);
			std::cout << _meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			boneSurface->DeepCopy(meshReader->GetOutput());	
		}
		
		// extracting outer surface
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		if (_meshFileNames[0].endsWith(".inp") || _meshFileNames[0].endsWith(".cdb")){
			surfaceExtractor->SetVolumeMesh(boneSurface);
			surfaceExtractor->Update();
			boneSurface->DeepCopy(surfaceExtractor->GetSurfaceMesh());	

		}
		
		/*
		RenderingMesh* renderingMesh2 = new RenderingMesh;
		renderingMesh2->SetMesh(boneSurface);
		renderingMesh2->SetRenderer(_ren);
		renderingMesh2->Update();
		
		
		// top
		RenderingOBB* obb = new RenderingOBB;
		obb->SetMesh(boneSurface);
		obb->CalculateOBB();
		obb->SetRenderer(_ren);
		obb->Update();
		obb->CalculateOBBtop();

		vtkPoints* top = vtkPoints::New();
		top = obb->GetPoints();
		RenderingOBB* obbHead = new RenderingOBB;
		obbHead->SetPoints(top);
		obbHead->CalculateOBB();
		obbHead->SetRenderer(_ren);
		obbHead->Update();
		std::cout << top->GetNumberOfPoints() << std::endl;

		// bottom
		RenderingOBB* obbBottom = new RenderingOBB;
		obbBottom->SetMesh(boneSurface);
		obbBottom->CalculateOBB();
		obbBottom->SetRenderer(_ren);
		obbBottom->Update();
		obbBottom->CalculateOBBbottom();

		vtkPoints* bottom = vtkPoints::New();
		bottom = obbBottom->GetPoints();
		RenderingOBB* renBottom = new RenderingOBB;
		renBottom->SetPoints(bottom);
		renBottom->CalculateOBB();
		renBottom->SetRenderer(_ren);
		renBottom->Update();
		std::cout << bottom->GetNumberOfPoints() << std::endl;
		*/
		
		
		// fitting inputs (mesh)
		MeshBoneImplantFitting* boneImplantFitting = new MeshBoneImplantFitting;
		boneImplantFitting->SetBoneMesh(boneSurface);
		boneImplantFitting->SetImplantMesh(implantReader->GetOutput());
		
		// implant-bone position
		if (_positionFlag == 1){
			
			std::cout << "proximal" << std::endl;
			boneImplantFitting->SetImplantPosition(_positionFlag);
			
		}
		else if (_positionFlag == 2){

			std::cout << "diaphyseal" << std::endl;
			boneImplantFitting->SetImplantPosition(_positionFlag);
		}
		else if (_positionFlag == 3){

			std::cout << "distal" << std::endl;
			boneImplantFitting->SetImplantPosition(_positionFlag);
		}
		else if ((_positionFlag != 1) && (_positionFlag != 2) && (_positionFlag != 3)){
			
			std::cout << "bone-implant position missing" << std::endl;
		}
	
		boneImplantFitting->ImplantPositionUpdate();


		// rendering
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(boneImplantFitting->GetImplantMesh());
		renderingMesh->SetRenderer(_ren);
		double color[3] = {221.0/255.0, 222.0/255.0, 223.0/255.0};
		renderingMesh->SetColor(color);
		renderingMesh->Update();
		
		RenderingOBB* implantOBB = new RenderingOBB;
		implantOBB->SetMesh(boneImplantFitting->GetImplantMesh());
		implantOBB->SetRenderer(_ren);
		implantOBB->CalculateOBB();
		implantOBB->Update();
		
		RenderingMesh* renderingMesh2 = new RenderingMesh;
		renderingMesh2->SetMesh(boneImplantFitting->GetBoneMesh());
		renderingMesh2->SetRenderer(_ren);
		renderingMesh2->Update();
		
		RenderingOBB* boneOBB = new RenderingOBB;
		boneOBB->SetMesh(boneImplantFitting->GetBoneMesh());
		boneOBB->SetRenderer(_ren);
		boneOBB->CalculateOBB();
		boneOBB->Update();
		

		
				
		boneImplantFitting->Update();
		
	}


}

/***************************** RENDERING ******************************/
void MainWindow::renderMetafile(){

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load image", " ", "*.mhd;;*.hdr");
	
	for (int i=0; i<fileNames.size(); i++){	
		std::cout << "rendering: " << fileNames[i].ascii() << std::endl;

		ImageHandler* imageHandler = new ImageHandler;
		imageHandler->SetImageFileName(fileNames[i]);
		imageHandler->MetafileReaderUpdate();
		imageHandler->SetThreshold(-800); // 
		imageHandler->MarchingCubesUpdate();
		
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(imageHandler->GetMesh());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();

		delete imageHandler;
		delete renderingMesh;
	}

}
void MainWindow::renderStl(){
	
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load mesh", " ", "*.stl");
	
	for (int i=0; i<fileNames.size(); i++){	
		std::cout << "rendering: " << fileNames[i].ascii() << std::endl;
			
		vtkSTLReader* stlReader = vtkSTLReader::New();
		stlReader->SetFileName(fileNames[i]);
		stlReader->Update();
		std::cout << "nodes: " << stlReader->GetOutput()->GetNumberOfPoints() << std::endl;
		std::cout << "elements: " <<  stlReader->GetOutput()->GetNumberOfCells() << std::endl;
					
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(stlReader->GetOutput());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();

		stlReader->Delete();
		delete renderingMesh;
	}
}
void MainWindow::renderSurfaceCdb(){

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load mesh", " ", "*.cdb");
	
	for (int i=0; i<fileNames.size(); i++){
		std::cout << "rendering: " << fileNames[i].ascii() << std::endl;
		
		MeshReaderMorpherSurface* reader = new MeshReaderMorpherSurface;
		reader->SetFileName(fileNames[i]);
		reader->Update();

		RenderingMesh* renderingMesh = new RenderingMesh;
		double color[3];
		color[0] = 232.0/255.0; color[1] = 184.0/255.0; color[2] = 45.0/255.0; // bone color
		//color[0] = 219.0/255.0; color[1] = 219.0/255.0; color[2] = 219.0/255.0;
		renderingMesh->SetColor(color);
		renderingMesh->SetMesh(reader->GetOutput());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();

		delete reader;
		delete renderingMesh;
	}

}
void MainWindow::renderVolumeCdb(){

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load mesh", " ", "*.cdb");
	
	for (int i=0; i<fileNames.size(); i++){
		std::cout << "rendering: " << fileNames[i].ascii() << std::endl;
		
		MeshReaderMorpherVolume* reader = new MeshReaderMorpherVolume;
		reader->SetFileName(fileNames[i]);
		reader->Update();
		
		std::cout << "extracting outer surface" << std::endl;
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		surfaceExtractor->SetVolumeMesh(reader->GetOutput());
		surfaceExtractor->Update();
		
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(surfaceExtractor->GetSurfaceMesh());
		double color[3];
		color[0] = 232.0/255.0; color[1] = 184.0/255.0; color[2] = 45.0/255.0; // bone color
		//color[0] = 219.0/255.0; color[1] = 219.0/255.0; color[2] = 219.0/255.0;
		//color[0] = 255.0/255.0; color[1] = 0.0/255.0; color[2] = 0.0/255.0;
		renderingMesh->SetColor(color);
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();

		delete reader;
		delete surfaceExtractor;
		delete renderingMesh;
	}

}
void MainWindow::renderVolumeInpAbaqus(){

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load mesh", " ", "*.inp");
	
	for (int i=0; i<fileNames.size(); i++){
		std::cout << "rendering: " << fileNames[i].ascii() << std::endl;
		
		MeshReaderAbaqus* reader = new MeshReaderAbaqus;
		reader->SetFileName(fileNames[i]);
		reader->Update();

		std::cout << "extracting outer surface" << std::endl;
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		surfaceExtractor->SetVolumeMesh(reader->GetOutput());
		surfaceExtractor->Update();
		
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(surfaceExtractor->GetSurfaceMesh());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();

		delete reader;
		delete surfaceExtractor;
		delete renderingMesh;
	}
	
	

}
void MainWindow::renderVtkPolyData(){
	
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load vtkPolyData", " ", "*.vtk");
	
	for (int i=0; i<fileNames.size(); i++){
		std::cout << "rendering: " << fileNames[i].ascii() << std::endl;
	
		vtkPolyDataReader* reader = vtkPolyDataReader::New();
		reader->SetFileName(fileNames[i].ascii());
		reader->Update();

		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(reader->GetOutput());
		renderingMesh->SetRenderer(_ren);
		double color[3]; color[0]=1.0; color[1]=0.0; color[2]=0.0;
		renderingMesh->SetColor(color);
		renderingMesh->Update();
	}
}
void MainWindow::renderPointsXyz(){

	QString fileName= QFileDialog::getOpenFileName(this, "Select files", " ", "*.txt"); 
	std::cout << "rendering: " << fileName.toAscii().data() << std::endl;
	
	PointReaderXyz* pointReader = new PointReaderXyz();
	pointReader->SetFileName(fileName);
	pointReader->Update();

	RenderingPoint* renderingPoint = new RenderingPoint;
	renderingPoint->SetPoints(pointReader->GetOutput());
	double color[3];
	//color[0] = 232.0/255.0; color[1] = 184.0/255.0; color[2] = 45.0/255.0; // bone color
	color[0] = 1.0; color[1] = 0.0; color[2] = 0.0;
	renderingPoint->SetColor(color);
	renderingPoint->SetRadius(5);
	renderingPoint->SetPointRatio(1);
	renderingPoint->SetRenderer(_ren);
	renderingPoint->Update();

	delete pointReader;
	delete renderingPoint;

}
void MainWindow::renderPointsXyzWithColorbar(){

	QString pointsFileName= QFileDialog::getOpenFileName(this, "Select point files", " ", "*.txt"); 
	QString colorbarFileName= QFileDialog::getOpenFileName(this, "Select colorbar files", " ", "*.txt"); 
	
	std::cout << "rendering points with colorbar" << std::endl;
	std::cout << pointsFileName.toAscii().data() << std::endl;
	std::cout << colorbarFileName.toAscii().data() << std::endl;

	PointReaderXyz* pointReader = new PointReaderXyz();
	pointReader->SetFileName(pointsFileName);
	pointReader->Update();

	VnlReaderVector* vectorReader = new VnlReaderVector;
	vectorReader->SetFileName(colorbarFileName);
	vectorReader->Update();
	vnl_vector<double> vnlColorbar;
	vnlColorbar = vectorReader->GetVnlVector();
	vtkDoubleArray* vtkColorbar = vtkDoubleArray::New();
	for (int i=0; i<vnlColorbar.size(); i++)
		vtkColorbar->InsertNextValue(vnlColorbar(i));

	RenderingPointWithColorbar* renderingPoint = new RenderingPointWithColorbar;
	renderingPoint->SetPoints(pointReader->GetOutput());
	renderingPoint->SetRadius(0.5);
	renderingPoint->SetPointRatio(1);
	renderingPoint->SetColorBar(vtkColorbar);
	//_background[0] = 198.0/255.0; _background[1] = 226.0/255.0; _background[2] = 255.0/255.0;
	_ren->SetBackground(_background); 
	renderingPoint->SetRenderer(_ren);
	renderingPoint->UpdateBW();

	delete pointReader;
	delete renderingPoint;
	delete vectorReader;

}
void MainWindow::renderPointsXyzWithColorbarColors(){

	QString pointsFileName= QFileDialog::getOpenFileName(this, "Select point files", " ", "*.txt"); 
	QString colorbarFileName= QFileDialog::getOpenFileName(this, "Select colorbar files", " ", "*.txt"); 
	
	std::cout << "rendering points with colorbar" << std::endl;
	std::cout << pointsFileName.toAscii().data() << std::endl;
	std::cout << colorbarFileName.toAscii().data() << std::endl;

	PointReaderXyz* pointReader = new PointReaderXyz();
	pointReader->SetFileName(pointsFileName);
	pointReader->Update();

	VnlReaderVector* vectorReader = new VnlReaderVector;
	vectorReader->SetFileName(colorbarFileName);
	vectorReader->Update();
	vnl_vector<double> vnlColorbar;
	vnlColorbar = vectorReader->GetVnlVector();
	vtkDoubleArray* vtkColorbar = vtkDoubleArray::New();
	for (int i=0; i<vnlColorbar.size(); i++)
		vtkColorbar->InsertNextValue(vnlColorbar(i));

	RenderingPointWithColorbar* renderingPoint = new RenderingPointWithColorbar;
	renderingPoint->SetPoints(pointReader->GetOutput());
	renderingPoint->SetRadius(0.5);
	renderingPoint->SetPointRatio(1);
	renderingPoint->SetColorBar(vtkColorbar);
	//_ren->SetBackground(_background); 
	renderingPoint->SetRenderer(_ren);
	renderingPoint->UpdateC();

	delete pointReader;
	delete renderingPoint;
	delete vectorReader;

}
void MainWindow::renderLine(){

	QString fileName= QFileDialog::getOpenFileName(this, "Select files", " ", "*.txt"); 
	std::cout << "rendering: " << fileName.toAscii().data() << std::endl;
	
	PointReaderXyz* pointReader = new PointReaderXyz();
	pointReader->SetFileName(fileName);
	pointReader->Update();

	vtkLineSource* line = vtkLineSource::New();
	double pt[3];
	pointReader->GetOutput()->GetPoint(0,pt);
	std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
	line->SetPoint1(pt);
	pointReader->GetOutput()->GetPoint(1,pt);
	std::cout << pt[0] << ' ' << pt[1] << ' ' << pt[2] << std::endl;
	line->SetPoint2(pt);
	vtkPolyDataMapper * mapper = vtkPolyDataMapper ::New();
	mapper->SetInput(line->GetOutput());
	vtkActor * actor = vtkActor ::New();
	actor->SetMapper(mapper);
	double color[3];
	color[0] = 1.0; color[1] = 0.0; color[2] = 0.0;
	actor->GetProperty()->SetColor(color);
	actor->GetProperty()->SetLineWidth(5.0);
	_ren->AddActor(actor);
	_ren->ResetCamera();
	_ren->GetRenderWindow()->Render();




}
void MainWindow::renderOBB(){

	QStringList meshFileNames = QFileDialog::getOpenFileNames(this, "Load bone meshes","C:/0.Data/test data/bone implant fitting","STL meshes(*.stl);;Abaqus quadratic volume meshes (*.inp);;Ansys quadratic volume meshes(*.cdb)");
	
	for (int i=0; i<meshFileNames.size(); i++){
		
		// read bone mesh
		vtkPolyData* boneSurface = vtkPolyData::New();
		
		// abaqus
		if (meshFileNames[i].endsWith(".inp")){
			MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
			meshReader->SetFileName(meshFileNames[i]);
			std::cout << meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			boneSurface->DeepCopy(meshReader->GetOutput());		
		}
		// ansys
		else if (meshFileNames[i].endsWith(".cdb")){
			MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
			meshReader->SetFileName(meshFileNames[i]);
			std::cout << meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			boneSurface->DeepCopy(meshReader->GetOutput());	
		}
		// stl
		else if (meshFileNames[i].endsWith(".stl")){
			vtkSTLReader* meshReader = vtkSTLReader::New();
			meshReader->SetFileName(meshFileNames[i]);
			std::cout << meshFileNames[i].ascii() << std::endl; 
			meshReader->Update();
			boneSurface->DeepCopy(meshReader->GetOutput());	
		}
		
		// extracting outer surface
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		if (meshFileNames[i].endsWith(".inp") || meshFileNames[i].endsWith(".cdb")){
			surfaceExtractor->SetVolumeMesh(boneSurface);
			surfaceExtractor->Update();
			boneSurface->DeepCopy(surfaceExtractor->GetSurfaceMesh());	

		}
		
		// rendering mesh
		RenderingMesh* renderingMesh2 = new RenderingMesh;
		renderingMesh2->SetMesh(boneSurface);
		renderingMesh2->SetRenderer(_ren);
		renderingMesh2->Update();
		
		// rendering OBB
		RenderingOBB* obb = new RenderingOBB;
		obb->SetMesh(boneSurface);
		obb->CalculateOBB();
		obb->SetRenderer(_ren);
		obb->Update();

		double one[3], two[3], three[3], four[3], five[3], six[3], seven[3], eight[3];
		obb->GetOBBverteces(one, two, three, four, five, six, seven, eight);
		std::cout << "Vertex numbering: considering a vertical parallelepiped" << std::endl;
		std::cout << "- bottom base, from the front left vertex, counter-clockwise: 1-3-7-4" << std::endl;
		std::cout << "- upper base, from the front left vertex, counter-clockwise: 2-6-8-5" << std::endl;
		std::cout << "1 (corner): " << one[0] << " " << one[1] << " " << one[2] << std::endl;
		std::cout << "2 (max): " << two[0] << " " << two[1] << " " << two[2] << std::endl;
		std::cout << "3 (mid): " << three[0] << " " << three[1] << " " << three[2] << std::endl;
		std::cout << "4 (min): " << four[0] << " " << four[1] << " " << four[2] << std::endl;
		double obbCenter[3];
		obb->GetOBBcenter(obbCenter);
		std::cout << "obb center: " << obbCenter[0] << " " << obbCenter[1] << " " << obbCenter[2] << std::endl;
		double centerOfMass[3];
		obb->GetMeshCenterOfMass(centerOfMass);
		std::cout << "center of mass: " << centerOfMass[0] << " " << centerOfMass[1] << " " << centerOfMass[2] << std::endl;
	
	}

}
void MainWindow::renderCoordinateSystem(){

	RenderingCoordinateSystem* renderingCS = new RenderingCoordinateSystem;
	renderingCS->SetRenderer(_ren);
	renderingCS->Update();

}
void MainWindow::euclideanDistancePoints(){

	QString fileName= QFileDialog::getOpenFileName(this, "Select the coordinate file", "C:/0.Data/test data/distance rendering/points", "*.txt"); 
		
	QStringList fileNames= QFileDialog::getOpenFileNames(this, "Select the two files", "C:/0.Data/test data/distance rendering/points", "*.txt"); 
	
	if (fileNames.size() > 2)
		std::cout << "only the first two inputs will be considered" << std::endl;

	std::cout << "rendering the difference between: " << std::endl;
	std::cout << fileNames[0].toAscii().data() << std::endl;
	std::cout << "and" << std::endl;
	std::cout << fileNames[1].toAscii().data() << std::endl;
	
	// coordinates
	PointReaderXyz* pointReader = new PointReaderXyz();
	pointReader->SetFileName(fileName);
	pointReader->Update();

	// differences for the toolbar
	VnlReaderVector* firstVectorReader = new VnlReaderVector;
	firstVectorReader->SetFileName(fileNames[0]);
	firstVectorReader->Update();

	VnlReaderVector* secondVectorReader = new VnlReaderVector;
	secondVectorReader->SetFileName(fileNames[1]);
	secondVectorReader->Update();

	StatisticsDistanceCalculator* distanceCalculator = new StatisticsDistanceCalculator;
	distanceCalculator->SetVectorOne(firstVectorReader->GetVnlVector());
	distanceCalculator->SetVectorTwo(secondVectorReader->GetVnlVector());
	distanceCalculator->CalculateEuclideanDistance();
	vnl_vector<double> vnlColorbar;
	vnlColorbar = distanceCalculator->GetDistanceVector();

	// rendering
	vtkDoubleArray* vtkColorbar = vtkDoubleArray::New();
	for (int i=0; i<vnlColorbar.size(); i++)
		vtkColorbar->InsertNextValue(vnlColorbar(i));

	RenderingPointWithColorbar* renderingPoint = new RenderingPointWithColorbar;
	renderingPoint->SetPoints(pointReader->GetOutput());
	renderingPoint->SetRadius(2);
	renderingPoint->SetColorBar(vtkColorbar);
	double colorRange[2]; colorRange[0] = 0.67; colorRange[1] = 0.0;
	renderingPoint->SetColorRange(colorRange);
	renderingPoint->SetRenderer(_ren);
	renderingPoint->Update();


}
void MainWindow::drawLineStop(){

	_cellPicker->RemoveAllObservers();
	
	// file name 
	if (_fileName.endsWith(".inp"))
		_fileName.replace(QString(".inp"), QString("_points.txt"));
	else if (_fileName.endsWith(".cdb"))
		_fileName.replace(QString(".cdb"), QString("_points.txt"));
	
	// print out
	for (int i=0; i<_coords->GetNumberOfPoints(); i++){
		double coord[3];
		_coords->GetPoint(i, coord);
		std::cout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << std::endl;
	}

	/*
	// save points
	PointWriterXyz* pointWriter = new PointWriterXyz;
	pointWriter->SetInput(_coords);
	pointWriter->SetFileName(_fileName.toAscii().data());
	pointWriter->Update();
	std::cout << ".txt load file written" << std::endl;
	*/
	if (_coords->GetNumberOfPoints() == 2){
		
		vtkLineSource* line = vtkLineSource::New();
		double pt[3];
		_coords->GetPoint(0,pt);
		line->SetPoint1(pt);
		_coords->GetPoint(1,pt);
		line->SetPoint2(pt);
		vtkPolyDataMapper * mapper = vtkPolyDataMapper ::New();
		mapper->SetInput(line->GetOutput());
		vtkActor * actor = vtkActor ::New();
		actor->SetMapper(mapper);
		double color[3];
		color[0] = 1.0; color[1] = 0.0; color[2] = 0.0;
		actor->GetProperty()->SetColor(color);
		actor->GetProperty()->SetLineWidth(5.0);
		_ren->AddActor(actor);
		
	
	}
	
	
}
void MainWindow::resetView(){

	_ren->RemoveAllViewProps();
	_ren->ResetCamera();
	_background[0] = 1.0; _background[1] = 1.0; _background[2] = 1.0;
	_ren->SetBackground(_background); // Background color white
  	qvtkWidget->GetRenderWindow()->Render();

	std::cout << "reset view" << std::endl;

}

/************************ FORMAT CONVERSION ***************************/
void MainWindow::dicomToMetafile(){

	std::cout << "image file conversion from .dcm to .mhd" << std::endl;
	std::cout << std::endl;

	QString folderName = QFileDialog::getExistingDirectory(this, "Select parent folder"); 
	
	QStringList list;
	QDir dir(folderName);
	QStringList filters;
	filters << "*SLICES";
	dir.setNameFilters(filters);
	list = dir.entryList();
		
	for(int i = 0; i < list.size(); ++i) {

		std::cout << (dir.absolutePath() + "/" + list[i]).toAscii().data() << std::endl;
		
		ImageHandler* imageHandler = new ImageHandler;
		
		imageHandler->SetDicomSerieFolderName(dir.absolutePath() + "/" + list[i]);
		imageHandler->DicomReaderUpdate();
				
		imageHandler->SetImageFileName(dir.absolutePath() + "/" + list[i] + ".mhd");
		imageHandler->MetafileWriterUpdate();
		
		delete imageHandler;
	}

	std::cout << "images converted" << std::endl;

}
void MainWindow::metafileToDicom(){

	std::cout << "image file conversion from .mhd to .dcm" << std::endl;
	std::cout << std::endl;

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load metafile image", " ", "*.mhd");
	
	for (int i=0; i<fileNames.size(); i++){

		std::cout << "Image: " << fileNames[i].ascii() << std::endl;
		
		ImageHandler* imageHandler = new ImageHandler;
		
		imageHandler->SetImageFileName(fileNames[i]);
		imageHandler->MetafileReaderUpdate();

		
		QDir dir(fileNames[i]);
		fileNames[i].replace(QString(".mhd"), QString(""));
		dir.mkdir(fileNames[i]);
		imageHandler->SetDicomSerieFolderName(fileNames[i]);
		imageHandler->DicomWriterUpdate();

		delete imageHandler;
	}

	std::cout << "images converted" << std::endl;
}
void MainWindow::metafileToAnalyse(){

	
	std::cout << "image file conversion from .mhd to .hdr" << std::endl;
	std::cout << std::endl;

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load metafile image", " ", "*.mhd");
	
	for (int i=0; i<fileNames.size(); i++){

		std::cout << "Image: " << fileNames[i].ascii() << std::endl;
		
		ImageHandler* imageHandler = new ImageHandler;
		
		imageHandler->SetImageFileName(fileNames[i]);
		imageHandler->MetafileReaderUpdate();
		
		fileNames[i].replace(QString(".mhd"), QString(".hdr"));	
		imageHandler->SetImageFileName(fileNames[i]);
		imageHandler->MetafileWriterUpdate();

		delete imageHandler;
	}

	std::cout << "images converted" << std::endl;
	
	/*
	typedef float VoxelType;
	static const unsigned int Dimension = 3;
	typedef itk::Image< VoxelType, Dimension > ImageType;
	typedef itk::ImageFileReader< ImageType > ImageReaderType;
	typedef itk::ImageFileWriter< ImageType > ImageWriterType;

	typedef float FieldVoxelType;
	typedef itk::Vector<FieldVoxelType,Dimension> VectorType;
	typedef itk::Image<VectorType,Dimension> FieldType; 
	typedef itk::ImageFileReader< FieldType > FieldReaderType;
	typedef itk::ImageFileWriter< FieldType > FieldWriterType;

	
	QString folderName = QFileDialog::getExistingDirectory(this, "Select folder"); 
	
	QStringList list;
	QDir dir(folderName);
	QStringList filters;
	

	// making new folders
	std::cout << "creating new folders" << std::endl;
	filters << "bone*"; 
	dir.setNameFilters(filters);
	list = dir.entryList();

	std::cout << "number of folders: " << list.size() << std::endl;
	
	for(int i = 0; i < list.size(); ++i) {
		
		QString folderName = (dir.absolutePath() + "/" + list[i]);
		std::cout << folderName.toAscii().data() << std::endl;
		QString boneNumber = list[i];
		boneNumber.remove("bone ");
		std::cout << boneNumber.toAscii().data() << std::endl;
			

		// instance
		QString imageName = folderName;
		imageName.append("\\instance.mhd");
		std::cout << "image name: " << imageName.toAscii().data() << std::endl;

		ImageReaderType::Pointer imageReader = ImageReaderType::New();
		imageReader->SetFileName(imageName.toAscii().data());
		imageReader->Update();

		imageName.replace(QString ("instance.mhd"), QString("instance_"));
		imageName.append(boneNumber);
		imageName.append(".mhd");
		std::cout << "new name: " << imageName.toAscii().data() << std::endl;
		
		ImageWriterType::Pointer imageWriter = ImageWriterType::New();
		imageWriter->SetFileName(imageName.toAscii().data());
		imageWriter->SetInput(imageReader->GetOutput());
		imageWriter->Update();

		// warped to reference
		imageName = folderName;
		imageName.append("\\int.mhd");
		std::cout << "image name: " << imageName.toAscii().data() << std::endl;

		imageReader->SetFileName(imageName.toAscii().data());
		imageReader->Update();

		imageName.replace(QString ("int.mhd"), QString("intensityInReference_"));
		imageName.append(boneNumber);
		imageName.append(".mhd");
		std::cout << "new name: " << imageName.toAscii().data() << std::endl;
		
		imageWriter->SetFileName(imageName.toAscii().data());
		imageWriter->SetInput(imageReader->GetOutput());
		imageWriter->Update();

		// dvf
		QString fieldName = folderName;
		fieldName.append("\\dvf.mhd");
		std::cout << "field name: " << fieldName.toAscii().data() << std::endl;

		FieldReaderType::Pointer fieldReader = FieldReaderType::New();
		fieldReader->SetFileName(fieldName.toAscii().data());
		fieldReader->Update();

		fieldName.replace(QString ("dvf.mhd"), QString("dvf_"));
		fieldName.append(boneNumber);
		fieldName.append(".mhd");
		std::cout << "new name: " << fieldName.toAscii().data() << std::endl;
		
		FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
		fieldWriter->SetFileName(fieldName.toAscii().data());
		fieldWriter->SetInput(fieldReader->GetOutput());
		fieldWriter->Update();
		
	}
	



	std::cout << "done" << std::endl;
	*/

	/*
	typedef signed short VoxelType;
	static const unsigned int Dimension = 3;
	typedef itk::Image< VoxelType, Dimension > ImageType;
	typedef itk::ImageFileReader< ImageHandler::ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageHandler::ImageType > WriterType;
	
	
	QString folderName = QFileDialog::getExistingDirectory(this, "Select folder"); 
	
	QStringList list;
	QDir dir(folderName);
	QStringList filters;
	

	// making new folders
	std::cout << "creating new folders" << std::endl;
	filters << "50*"; 
	dir.setNameFilters(filters);
	list = dir.entryList();

	std::cout << "number of folders: " << list.size() << std::endl;
	
	for(int i = 0; i < list.size(); ++i) {
		
		QString folderName = (dir.absolutePath() + "/" + list[i]);
		std::cout << folderName.toAscii().data() << std::endl;

		folderName.append(QString("-xx-x-xxx-xxx-xxx_mandible-l-Segmented-SLICES"));
		std::cout << folderName.toAscii().data() << std::endl;
		dir.mkdir(folderName);
		
		folderName.replace(QString("Segmented"), QString("Classified"));
		std::cout << folderName.toAscii().data() << std::endl;
		dir.mkdir(folderName);
	}
	

	// filling in the new folders
	std::cout << "filling in the new folders" << std::endl;
	
	for(int i = 0; i <list.size(); ++i) {

		ReaderType::Pointer reader = ReaderType::New();
		WriterType::Pointer writer = WriterType::New();
			
		QString folderName = (dir.absolutePath() + "/" + list[i]);
		
		// Segmented folders
		QString fileName = dir.absolutePath() + "/" + list[i] + "/" + "Anonymized" + list[i] + ".hdr";
		std::cout << fileName.toAscii().data() << std::endl;
		reader->SetFileName(fileName.toAscii().data());
		reader->Update();
		
		QString destinationFileName = folderName;
		destinationFileName.append(QString("-xx-x-xxx-xxx-xxx_mandible-l-Segmented-SLICES"));
		destinationFileName.append("/");
		destinationFileName.append("image.hdr");
		std::cout << destinationFileName.toAscii().data() << std::endl;
		writer->SetFileName(destinationFileName.toAscii().data());
		writer->SetInput(reader->GetOutput());
		writer->Update();

		// Classified folders
		fileName.replace(QString(".hdr"), QString(".Labels.hdr"));
		std::cout << fileName.toAscii().data() << std::endl;
		reader->SetFileName(fileName.toAscii().data());
		reader->Update();

		destinationFileName.replace(QString("Segmented"), QString("Classified"));
		std::cout << destinationFileName.toAscii().data() << std::endl;
		writer->SetFileName(destinationFileName.toAscii().data());
		writer->SetInput(reader->GetOutput());
		writer->Update();
	}

	std::cout << "done" << std::endl;
	*/

}
void MainWindow::morpherSurfaceCdbToStl(){

	std::cout << "surface mesh file conversion from .cdb to .stl" << std::endl;
	std::cout << std::endl;
	
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load Ansys Morpher surface meshes", " ", "*.cdb");
	 
	for (int i=0; i<fileNames.size(); i++){
	
		std::cout << "Mesh: " << fileNames[i].ascii() << std::endl;
		
		MeshReaderMorpherSurface* meshReaderMorpherSurface = new MeshReaderMorpherSurface;
		meshReaderMorpherSurface->SetFileName(fileNames[i]);
		meshReaderMorpherSurface->Update();

		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(meshReaderMorpherSurface->GetOutput());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();

		vtkSTLWriter* stlWriter = vtkSTLWriter::New();
		fileNames[i].replace(QString(".cdb"), QString(".stl"));
		stlWriter->SetFileName(fileNames[i].ascii());
		stlWriter->SetInput(meshReaderMorpherSurface->GetOutput());
		stlWriter->Update();

		delete meshReaderMorpherSurface;
		delete renderingMesh;
		stlWriter->Delete();

	}
		
}
void MainWindow::morpherVolumeCdbToVtk(){

	std::cout << "volume mesh file conversion from .cdb to .vtk" << std::endl;
	std::cout << std::endl;

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load Ansys Morpher volume meshes", " ", "*.cdb");
	 
	for (int i=0; i<fileNames.size(); i++){
	
		std::cout << "Mesh: " << fileNames[i].ascii() << std::endl;
		
		MeshReaderMorpherVolume* meshReaderMorpherVolume = new MeshReaderMorpherVolume;
		meshReaderMorpherVolume->SetFileName(fileNames[i]);
		meshReaderMorpherVolume->Update();

		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(meshReaderMorpherVolume->GetOutput());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();
		
		vtkPolyDataWriter* polyDataWriter = vtkPolyDataWriter::New();
		fileNames[i].replace(QString(".cdb"), QString(".vtk"));
		polyDataWriter->SetFileName(fileNames[i].ascii());
		polyDataWriter->SetInput(meshReaderMorpherVolume->GetOutput());
		polyDataWriter->Update();

		delete meshReaderMorpherVolume;
		delete renderingMesh;
		polyDataWriter->Delete();
	}
		
}
void MainWindow::abaqusToCdb(){

	
	std::cout << "volume mesh file conversion from .inp to .cdb" << std::endl;
	std::cout << std::endl;

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load abaqus format mesh file", " ", "*.inp");
	 
	for (int i=0; i<fileNames.size(); i++){
	
		std::cout << "Mesh: " << fileNames[i].ascii() << std::endl;

		MeshReaderAbaqus* reader = new MeshReaderAbaqus;
		reader->SetFileName(fileNames[i]);
		reader->Update();

		MeshWriterAnsys* writer = new MeshWriterAnsys;
		fileNames[i].replace(QString(".inp"), QString(".cdb"));
		writer->SetFileName(fileNames[i]);
		writer->MeshOn();
		writer->SetMesh(reader->GetOutput());
		writer->Update();
		std::cout << ".cbd file written" << std::endl;
		
		std::cout << "extracting outer surface for rendering" << std::endl;
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		surfaceExtractor->SetVolumeMesh(reader->GetOutput());
		surfaceExtractor->Update();
		
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(surfaceExtractor->GetSurfaceMesh());
		renderingMesh->SetRenderer(_ren);
		//renderingMesh->Update();
		
		std::cout << "file converted" << std::endl;

		delete reader;
		delete writer;
		delete surfaceExtractor;
		delete renderingMesh;

		
	}
}
void MainWindow::cdbToInpAbaqus(){

	/*
	// neutral format (netgen to abaqus)
	std::cout << "volume mesh file conversion from .cdb to .inp(abaqus)" << std::endl;
	std::cout << std::endl;

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load neutral format volume meshes", " ", "*.inp");
		 
	for (int i=0; i<fileNames.size(); i++){
	
		std::cout << "Mesh: " << fileNames[i].ascii() << std::endl;

		MeshReaderNeutralFormat* reader = new MeshReaderNeutralFormat;
		reader->SetFileName(fileNames[i]);
		reader->Update();

		MeshWriterAbaqus* writer = new MeshWriterAbaqus;
		fileNames[i].replace(QString(".inp"), QString("2.inp"));
		writer->SetFileName(fileNames[i]);
		writer->MeshOn();
		writer->SetMesh(reader->GetOutput());
		writer->Update();
		std::cout << ".inp file written" << std::endl;
		
		std::cout << "extracting outer surface for rendering" << std::endl;
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		surfaceExtractor->SetVolumeMesh(reader->GetOutput());
		surfaceExtractor->Update();
		
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(surfaceExtractor->GetSurfaceMesh());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();
		
		std::cout << "file converted" << std::endl;

		delete reader;
		delete writer;
		delete surfaceExtractor;
		delete renderingMesh;
	}
	*/


	std::cout << "volume mesh file conversion from .cdb to .inp(abaqus)" << std::endl;
	std::cout << std::endl;

	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load Ansys Morpher volume meshes", " ", "*.cdb");
		 
	for (int i=0; i<fileNames.size(); i++){
	
		std::cout << "Mesh: " << fileNames[i].ascii() << std::endl;

		MeshReaderMorpherVolume* reader = new MeshReaderMorpherVolume;
		reader->SetFileName(fileNames[i]);
		reader->Update();

		MeshWriterAbaqus* writer = new MeshWriterAbaqus;
		fileNames[i].replace(QString(".cdb"), QString(".inp"));
		writer->SetFileName(fileNames[i]);
		writer->MeshOn();
		writer->SetMesh(reader->GetOutput());
		writer->Update();
		std::cout << ".inp file written" << std::endl;
		
		std::cout << "extracting outer surface for rendering" << std::endl;
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		surfaceExtractor->SetVolumeMesh(reader->GetOutput());
		surfaceExtractor->Update();
		
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(surfaceExtractor->GetSurfaceMesh());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();
		
		std::cout << "file converted" << std::endl;

		delete reader;
		delete writer;
		delete surfaceExtractor;
		delete renderingMesh;
	}
	
}
void MainWindow::ptStlToTxt(){

	std::cout << "point file conversion from .stl to .txt" << std::endl;
	std::cout << std::endl;
	
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load .stl files", " ", "*.stl");

	for (int i=0; i<fileNames.size(); i++){
		
		vtkSTLReader* stlReader = vtkSTLReader::New();
		stlReader->SetFileName(fileNames[i].ascii());
		stlReader->Update();
		std::cout << "nodes: " << stlReader->GetOutput()->GetNumberOfPoints() << std::endl;
		std::cout << "elements: " <<  stlReader->GetOutput()->GetNumberOfCells() << std::endl;
			

		fileNames[i].replace(QString(".stl"), QString("_nodes.txt"));

		PointWriterXyz* pointWriter = new PointWriterXyz;
		pointWriter->SetInput(stlReader->GetOutput()->GetPoints());
		pointWriter->SetFileName(fileNames[i]);
		pointWriter->Update();

		stlReader->Delete();
		delete pointWriter;
	}
}
void MainWindow::ptVtkToTxt(){

	std::cout << "point file conversion from .vtk to .txt" << std::endl;
	std::cout << std::endl;
	
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Load .vtk files", " ", "*.vtk");

	for (int i=0; i<fileNames.size(); i++){
		
		vtkPolyDataReader* reader = vtkPolyDataReader::New();
		reader->SetFileName(fileNames[i].ascii());
		reader->Update();
		
		fileNames[i].replace(QString(".vtk"), QString("_nodes.txt"));
	
		PointWriterXyz* pointWriter = new PointWriterXyz;
		pointWriter->SetInput(reader->GetOutput()->GetPoints());
		pointWriter->SetFileName(fileNames[i]);
		pointWriter->Update();

		reader->Delete();
		delete pointWriter;
	}
}

/*
void MainWindow::propagationForMorpher(){

	//this function is not connected to Tools since the propagation of the reference landmarsk 
	// to the moving bones is not needed
	QStringList fileNames = QFileDialog::getOpenFileNames(this, "Select files", "/home", "*.txt *velocity_field.mhd"); 
	
	MeshMorpherPointPropagation* morpher = new MeshMorpherPointPropagation;
	morpher->SetFileNames(fileNames);
	morpher->setRenderer(_ren);
	morpher->Update();
}
*/


/*
void MainWindow::dataProjection(){

	_fileNames = QFileDialog::getOpenFileNames(this, "Load volume meshes", " ", "Ansys volume mesh (*.cdb)");
	QStringList eigenvectorFileNames = QFileDialog::getOpenFileNames(this, "Load eigenvectors", " ", "(*.vtk)");
	
	PCA3dProjectionMesh* projection = new PCA3dProjectionMesh;
	projection->SetFileNames(_fileNames);
	projection->SetEigenvectorsFileNames(eigenvectorFileNames);
	projection->Update();

	RenderingPoint* rendering = new RenderingPoint;
	rendering->SetPoints(projection->GetProjectedCoords());
	rendering->SetRadius(20);
	rendering->SetRenderer(_ren);
	rendering->Update();


}
*/

/*
void MainWindow::loadMeshBC(){

	_fileName = QFileDialog::getOpenFileName(this, "Load bone meshes","C:/0.Data/1. Registered data/0. ref 00705","Abaqus quadratic volume meshes (*.inp);;Ansys quadratic volume meshes(*.cdb)");
	_boneSurface = vtkPolyData::New();
		
		// abaqus
		if (_fileName.endsWith(".inp")){
			MeshReaderAbaqus* meshReader = new MeshReaderAbaqus;
			meshReader->SetFileName(_fileName);
			std::cout << _fileName.toAscii().data() << std::endl; 
			meshReader->Update();
			_boneSurface->DeepCopy(meshReader->GetOutput());		
		}
		// ansys
		else if (_fileName.endsWith(".cdb")){
			MeshReaderMorpherVolume* meshReader = new MeshReaderMorpherVolume;
			meshReader->SetFileName(_fileName);
			std::cout << _fileName.toAscii().data() << std::endl; 
			meshReader->Update();
			_boneSurface->DeepCopy(meshReader->GetOutput());	
		}

		// extracting outer surface
		MeshExtractOuterSurface* surfaceExtractor = new MeshExtractOuterSurface;
		surfaceExtractor->SetVolumeMesh(_boneSurface);
		surfaceExtractor->Update();
		
		RenderingMesh* renderingMesh = new RenderingMesh;
		renderingMesh->SetMesh(surfaceExtractor->GetSurfaceMesh());
		renderingMesh->SetRenderer(_ren);
		renderingMesh->Update();

}

void MainWindow::magnitudeL(double value){

	_loadMagnitude = value;
	std::cout << "Load magnitude: " << _loadMagnitude << " N" <<  std::endl;
	
}
void MainWindow::pickLPoints(){

	// picker
	std::cout << "pick first the load application point and then the direction one" << std::endl;
	std::cout << "(position the mouse on the point + press p)" << std::endl;
	_pointPicker->SetMesh(_boneSurface);
	_pointPicker->SetRenderer(_ren);
	double color[3];
	color[0] = 0.0; color[1] = 0.0; color[2] = 0.8;
	_pointPicker->SetColor(color);
	_pointPicker->SetPicker(_cellPicker);
	_coords = vtkPoints::New();
	_pointPicker->SetPoints(_coords);
	_pointPicker->Update();
	
}

void MainWindow::saveL(){

	_cellPicker->RemoveAllObservers();
	
	// force direction line drawing
	if (_coords->GetNumberOfPoints() == 2){
		
		vtkLineSource* line = vtkLineSource::New();
		double pt[3];
		_coords->GetPoint(0,pt);
		line->SetPoint1(pt);
		_coords->GetPoint(1,pt);
		line->SetPoint2(pt);
		vtkPolyDataMapper * mapper = vtkPolyDataMapper ::New();
		mapper->SetInput(line->GetOutput());
		vtkActor * actor = vtkActor ::New();
		actor->SetMapper(mapper);
		double color[3];
		color[0] = 0.0; color[1] = 0.0; color[2] = 0.8;
		actor->GetProperty()->SetColor(color);
		actor->GetProperty()->SetLineWidth(5.0);
		_ren->AddActor(actor);
	
	}
	
	// coordinates
	_coords = _pointPicker->GetPoints();
	for (int i=0; i<_coords->GetNumberOfPoints(); i++){
		double coord[3];
		_coords->GetPoint(i, coord);
	}

	// load
	vtkDoubleArray* loadMagnitude = vtkDoubleArray::New();
	loadMagnitude->SetNumberOfValues(2);
	loadMagnitude->InsertValue(0, _loadMagnitude);
	loadMagnitude->InsertValue(1, _loadMagnitude);

	// saving 
	PointWriterXyzId* pointWriter = new PointWriterXyzId; // point writer ID
	pointWriter->SetInput(_coords);
	pointWriter->SetIdVector(loadMagnitude);

	if (_fileName.endsWith(".inp"))
		_fileName.replace(QString(".inp"), QString("_L.txt"));
	else if (_fileName.endsWith(".cdb"))
		_fileName.replace(QString(".cdb"), QString("_L.txt"));
	pointWriter->SetFileName(_fileName.toAscii().data());
	pointWriter->Update();
	std::cout << ".txt load file written" << std::endl;

}
void MainWindow::typeBC(){

	_BCtype = 1; // encastre

}
void MainWindow::pickBCPoints(){

	// picker
	std::cout << "pick BC application points (position the mouse on the point + press p)" << std::endl;
	_pointPicker->SetMesh(_boneSurface);
	_pointPicker->SetRenderer(_ren);
	double color[3];
	color[0] = 0.0; color[1] = 0.8; color[2] = 0.0;
	_pointPicker->SetColor(color);
	_pointPicker->SetPicker(_cellPicker);
	_coords = vtkPoints::New();
	_pointPicker->SetPoints(_coords);
	_pointPicker->Update();

}
void MainWindow::saveBC(){

	_cellPicker->RemoveAllObservers();
	
	// coordinates
	_coords = _pointPicker->GetPoints();
	for (int i=0; i<_coords->GetNumberOfPoints(); i++){
		double coord[3];
		_coords->GetPoint(i, coord);
	}

	// BC type
	vtkDoubleArray* typeBC = vtkDoubleArray::New();
	typeBC->SetNumberOfValues(_coords->GetNumberOfPoints());
	for (int i=0; i<_coords->GetNumberOfPoints();i++)
		typeBC->InsertValue(i, _BCtype);
	
	// saving 
	PointWriterXyzId* pointWriter = new PointWriterXyzId; // point writer ID
	pointWriter->SetInput(_coords);
	pointWriter->SetIdVector(typeBC);

	if (_fileName.endsWith(".inp"))
		_fileName.replace(QString(".inp"), QString("_BC.txt"));
	else if (_fileName.endsWith(".cdb"))
		_fileName.replace(QString(".cdb"), QString("_BC.txt"));
	pointWriter->SetFileName(_fileName.toAscii().data());
	pointWriter->Update();
	std::cout << ".txt boundary condition file written" << std::endl;
}
*/