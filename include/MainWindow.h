/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <Ui_MainWindow.h>

#include <RenderingPointPicker.h>

#include <itkImage.h>

#include <vtkCellPicker.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkPointPicker.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include <QListWidget>
#include <QStringList>

class MainWindow;
class QVTKInteractor;

using namespace rendering;


class DataPacked
{
public:
	MainWindow* mainWindow;
};

/*
* Class linked to the GUI "tools" created with Qt
*/

class MainWindow : public QMainWindow, private Ui_MainWindow
{
    Q_OBJECT

public:
   
	// Constructor/Destructor
    MainWindow(QWidget* parent = 0);
    ~MainWindow();
	
	// scalar image
	typedef signed short VoxelType;
	static const unsigned int Dimension = 3;
	typedef itk::Image< VoxelType, Dimension > ImageType;
	ImageType::Pointer _image;
	
	// field image
	typedef float FieldVoxelType;
	typedef itk::Vector<FieldVoxelType,Dimension> VectorType;
	typedef itk::Image<VectorType,Dimension> FieldType; 
	FieldType::Pointer _field;

	// mesh
	vtkPolyData* _mesh;

	//itk::CenteredAffineTransform< double,3 > AffineTransformType;
	


public slots:
	
	/************************* Image tool - slots *************************/ 
	
	// ad hoc preprocessing

	/**
	* Loads a list of .mhd images and resamples them
	* ! Set the new spacing in ImageHandler constructor  
	*/
	void resample();
	/**
	* Loads one .mhd image and transforms it to vtk coordinate systems
	*/
	void itkToVtkToItk();
	/**
	* Executes a pseudo calibration on CT images. 0 HU are put to 0.9937mg/mm3, the max HU is put to 1.5mg/mm3
	*/
	void pseudoCalibration();
	/**
	* Loads a list of .mhd masked images which have three label layers (100, 200, 300) with transitions and eliminates the transitions
	*/
	void maskGreyLevel();
	/**
	* Loads one .mhd image and its .mhd correspondent mask and writes the .mhd masked image 
	* refere to MaskToImage() in ImageHandler.h (when the mask is 0, the original image is put to -1024 or to the image background value)
	*/
	void maskImage();
	/**
	* Loads folders of .dcm images and masks and saves the masked images as .dcm
	* refere to MaskToImage() in ImageHandler.h (when the mask is 0, the original image is put to -1024 or to the image background value)
	* ! Loads the parent folder (i.e. the folder that contains all the .dcm image and mask folders)
	* ! Set the folder file name filter manually
	* ! Create in the parent folder a folder called "ReSegmented" that will contain all the masked .dcm image folders
	*/
	void binaryMaskDcm();

	
	// correspondeces - registration

	/**
	* Finds the vf that is closest to the average one. See SelectReferenceAsMinimumDistanceToAverage() from VectorImageHandler.h
	*/
	void findReference();
	/**
	* Loads a list of .mhd velocity fields and transforms them in deformation vector fields
	* It uses VFtoDVF() from VectorImageHandler.h 
	*/
	void vfToDvf();
	/**
	* Loads a list of .mhd velocity fields and the .mhd reference image. VF are transformed to DVF and then applied to the reference in ordert to warp it to the moving shape
	*/
	void recreateMoving();
	/**
	* Loads a list of .mhd velocity fields and transforms them in deformation vector fields
	* It uses VFtoDVFInverted() from VectorImageHandler.h 
	*/
	void vfToDvfInverted();
	/**
	* Loads a list of .mhd velocity fields and the correspondent .mhd images. 
	* DVFs are applied to the correspondent moving in order to obtain images that have the shape of the reference and the intensity of the moving.
	* The output images are used for PCAImageCalculation().
	*/
	void recreateReference();
	
	
	
	/**
	* 
	*/
	void surfaceComparisonImage();
	/**
	* 
	*/
	void volumeComparisonImage();
	
	
	// statistical model - calculation
	/**
	* Flags if the shape model for images has to be calculated. It is a deformation model. 
	*/
	void shapeImagePCA();
	/**
	* Flags if the intensity model for images has to be calculated. It is an intensity model.
	*/
	void intensityImagePCA();
	/**
	* Flags if the combined model for images has to be calculated. 
	*/
	void combinedImagePCA();
	/**
	* Calculates PCA, based on the flag previously choosen.
	*/
	void calculateImagePCA();

	// statistical model - validation
	/**
	* Recreates instances (for testing model representation and generation abilities)
	*/
	void recreateImageInstance();
	
	

	// instance creation
	/**
	* Number of images to be created
	*/
	void nOfImageInstances(int value);
	/**
	* Number of modes to use to create instances
	*/
	void nOfImageModes(int value);
	/**
	* Creates images
	*/
	void createImageInstances();
	


	/************************* Mesh tool - slots **************************/
	
	// creation
	
	/**
	* Loads a list of .mhd or .stl files
	*/
	void loadImageMesh();
	/**
	* Gets the simplify parameter
	*/
	void simplify(double value);
	/**
	* Gets the smooth number of iterations
	*/
	void smooth(int value);
	/**
	* Creates meshes both form images (.mhd) and meshes (.stl), simplifying, smoothing and running MRFSurface
	* ! In case of images, set the wanted marching cubes threshold manually
	*/
	void createMesh();
	/**
	* Loads a list of .stl meshes
	*/
	void loadSurfaceMesh();
	/**
	* Creates volume meshes using NetGen
	* ! Set the wanted parameters manually
	*/
	void createVolumeMesh();
	/**
	* Evaluates the jacobian of the meshes. Prints out if there are some 0 or negative volumes.
	*/
	void meshJacobian();
	/**
	* Evaluates the edge ratio of the meshes.
	*/
	void meshEdgeRatio();
	/**
	* Evaluates the angle of the meshes. 
	*/
	void meshAngle();


	// morpher processing

	/**
	* Loads one .stl mesh, renders it and allows to pick the wanted landmarks
	*/
	void pickLandmarks();
	/**
	* Loads the reference x y z ID .txt files and writes the x y z .txt file of the picked landmarks and the ansys morpher .xml file
	*/
	void createXml();
	/**
	* Loads the previously created .xml file and the new reference x y z ID .txt files and overwrites the new reference IDs
	*/
	void createXmlFromFiles();
	/**
	* Loads a list of .xml files and writes the ansys morpher
	*/
	void runMorpher();
	/**
	* Executes the surface validation. See RegistrationValidationSurfaceMesh.h for details.
	*/
	void surfaceValidation();
	/**
	* Executes the volume validation. See RegistrationValidationVolumeMesh.h for details.
	*/
	void volumeValidation();

	
	// statistical model
	/**
	* Extracts nodes coordinates and grey levels from .cdb files (to be extended) and writes them as vectors
	*/
	void extractCoordinatesAndIntensities();
	/**
	* Extracts nodes coordinates and elements
	*/
	void extractNodesAndElements();
	/**
	* Aligns the mesh using the RegistrationProcrustesAlignment.h class. Only .cdb files (to be extended)
	*/
	void meshAlignment();
	/**
	* Selects the shape pca
	*/
	void shapeMeshPCA();
	/**
	* Selects the intensity pca
	*/
	void instensityMeshPCA();
	/**
	* Selects the combined pca
	*/
	void combinedMeshPCA();
	/**
	* Executes the chosen pca
	*/
	void calculateMeshPCA();
	/**
	* Recreates instances (for testing model representation and generation abilities)
	*/
	void recreateMeshInstance();
	

	// instance creation
	/**
	* Number of meshes to be created
	*/
	void nOfMeshInstances(int value);
	/**
	* Number of modes to be used
	*/
	void nOfMeshModes(int value);
	/**
	* Creates meshes
	*/
	void createMeshInstances();



	/***************************** FEM tool ******************************/ 
	
	// material properties
	
	/**
	* Loads the Abaqus quadratic tetra mesh
	*/
	void loadMesh();
	/**
	* Loads the .mhd image
	*/
	void loadImage();
	/**
	* Assigns the first term of the law
	*/
	void assignmentLawFirst(double value);
	/**
	* Assigns the second term of the law
	*/
	void assignmentLawSecond(double value);
	/**
	* Assignment to mesh elements
	*/
	void assignToElements();
	/**
	* Assignment to mesh nodes
	*/
	void assignToNodes();
	/**
	* Creation of an abaqus output file
	*/
	void abaqusMatProp();
	/**
	* Creation of an ansys output file (to do)
	*/
	void txtMatProp();
	/**
	* Assing mechanical properties and save the file
	*/
	void assign();

	
	// boundary conditions
	
	/**
	* Picks the force and boundary condition application points
	*/
	void pickApplicationPoints();
	/**
	* Saves the picked application points
	*/
	void saveApplicationPoints();
	/**
	* Propagates the application points
	*/
	void propagateApplicationPoints();
	/**
	* Calculates femur lenghts, body height, BMI, body weight, forces, force and BC coordinate system unit vectors
	*/
	void walking();
	/**
	* Calculates forces and boundary condition for the sideways fall
	*/
	void falling();
	/**
	* 
	*/
	void extractNeckIntensities();
	


	// finite element simulations
	
	/**
	* Creates the abaqus file for walking
	*/
	void createAbaqusInp();
	/**
	* Creates the abaqus file for walking
	*/
	void createAbaqusInpFalling();
	
	/**
	* Runs the abaqus simulation
	*/
	void runSimulation();

	
	// implant fitting
	
	/**
	* Loads the bone .inp or .cdb
	*/
	void loadBoneMesh();
	/**
	* Loads the implant .stl mesh
	*/
	void loadImplantMesh();
	/**
	* Sets bone-implant proximal position 
	*/
	void proximal();
	/**
	* Sets bone-implant diaphyseal position
	*/
	void diaphyseal();
	/**
	* Sets bone-implant distal position
	*/
	void distal();
	/**
	* Computes the fitting and saves new mesh position and distances 
	*/
	void computeAndSave();
	


	/***************************** Rendering ******************************/
	/**
	* Loads a list of metafiles, extracts the mesh using MarchingCubesUpdate() and renders it
	* ! Set the wanted marching cubes threshold manually
	*/
	void renderMetafile();
	/**
	* Loads a list of .stl files and renders them
	*/
	void renderStl();
	/**
	* Loads a list of ansys .cdb files (quadratic triangles) and renders them
	*/
	void renderSurfaceCdb();
	/**
	* Loads a list of ansys .cdb files (quadratic tetras), extracts the outer surface (linear triangles) and renders them
	*/
	void renderVolumeCdb();
	/**
	* Loads a list of abaqus .inp files (quadratic tetras), extracts the outer surface (linear triangles) and renders them
	*/
	void renderVolumeInpAbaqus();
	/**
	* Loads a list of vtkPolyData (*.vtk) and renders them
	*/
	void renderVtkPolyData();
	/**
	* Loads a .txt file of x y z points and renders them as spheres
	* ! Set the wanted sphere radius manually
	*/
	void renderPointsXyz();
	/**
	* Loads a .txt file of x y z points and a .txt file of a double vector and renders them as black and white spheres 
	* ! Set the wanted sphere radius manually
	*/
	void renderPointsXyzWithColorbar();	
	/**
	* Loads a .txt file of x y z points and a .txt file of a double vector and renders them as coloured spheres 
	* ! Set the wanted sphere radius manually
	*/
	void renderPointsXyzWithColorbarColors();	
	/**
	* Renders a line given 2 pt
	*/
	void renderLine();
	/**
	* Renders the OBB of the loaded mesh
	*/
	void renderOBB();
	/**
	* Renders the coordinate system in 0,0,0
	*/
	void renderCoordinateSystem();
	/**
	* 
	*/
	void drawLineStop();
	
	
	/**
	* Renders the euclidean difference of two pointsets with a color toolbar
	*/
	void euclideanDistancePoints();
	/**
	* Resets the window view, not the camera position
	*/
	void resetView();

	// Format conversion
	/**
	* Loads folders of .dcm images and tranforms them in .mhd files 
	* ! Loads the parent folder (i.e. the folder that contains all the .dcm image folders to be converted)
	* ! Set the folder file name filter manually
	*/
	void dicomToMetafile();
	/**
	* Loads a list of .mhd files and transforms them in dicom folder files
	*/
	void metafileToDicom();
	/**
	* Loads a list of .mhd files and tranforms them in .hdr
	*/
	void metafileToAnalyse();
	/**
	* Loads a list of surface ansys .cdb files (ansys morpher format) and transforms them in .stl 
	*/
	void morpherSurfaceCdbToStl();
	/**
	* Loads a list of volume ansys .cdb files (ansys morpher format) and transforms them in .vtk polydata 
	*/
	void morpherVolumeCdbToVtk();
	/**
	* Loads a list of volume quadratic tetra abaqus .inp files and transforms them in .cdb
	*/
	void abaqusToCdb();
	/**
	* Loads a list of volume quadratic tetra ansys .db files and transforms them in abaqus .inp
	*/
	void cdbToInpAbaqus();
	/**
	* Loads a list of .stl meshes, extracts the nodes and writes them in a x y z .txt file
	*/
	void ptStlToTxt();
	/**
	* Loads a list of .vtk meshes, extracts the nodes and writes them in a x y z .txt file
	*/
	void ptVtkToTxt();
	
	
private:

	// vtkWidget
	vtkRenderer* _ren;
	double _background[3];
	double _simplifyValue;
	int _smoothNofIteration;

	// fileNames string
	QStringList _fileNames;
	QString _fileName;
	QStringList _meshFileNames;
	QStringList _imageFileNames;

	// Data member - Image tool
	int _nOfImageInstances;
	int _nOfImageModes;

	// Data member - Mesh tool
	RenderingPointPicker* _pointPicker; // landmark picker
	vtkCellPicker* _cellPicker; 
	vtkPoints* _coords;
	double _modelVariability;
	int _nOfInstances;
	int _modeInterval;
	bool _flagMeshEigenvalue;
	bool _flagMeshEigenvector;
	bool _flagMeshAverage;
	double _PCAflag;
	int _nOfMeshInstances;
	int _nOfMeshModes;

	// Data member - FEM tool
	double _assignmentLawOne;
	double _assignmentLawTwo;
	int _nodeElementFlag;
	int _fileType;
	vtkPolyData* _boneSurface;
	int _positionFlag;	
	double _loadMagnitude;
	int _BCtype;
	bool _meshFlag;
	bool _mechPropFlag;
	bool _loadFlag;
	bool _bCInp;


	
};

#endif // MAINWINDOW_H