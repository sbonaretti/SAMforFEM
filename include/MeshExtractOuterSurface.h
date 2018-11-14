/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#ifndef	MESHEXTRACTOUTERSURFACE_H
#define MESHEXTRACTOUTERSURFACE_H

#include <vtkPolyData.h>

#include <vector>
#include <map>

namespace mesh{


	/**
	* It extracts the outer surface from a volume mesh. It extract only the vertices of the triangles (not the intermediate points in case of quadratic tetras).
	* The algorithm is based on the criterion that the surface faces are not shared between elements.
	*
	* vtkGeometryFilter and the filters of the same family do not work properly neither when considering the quadratic tetra mesh
	* as vtkPolyData (cell_type=7, i.e. vtk polygon) nor when as vtkUnstructuredGrid (cell_type=24, i.e. vtk_quadratic_tetra)
	*
	* The outer surface faces don not have the normals with the same direction.
	*
	* The class adapts the code of FEMAssigner by Andreas Siegrist to vtkPolyData.
	*/

	class MeshExtractOuterSurface{
		
	public:
		
		// constructor/destructor
		MeshExtractOuterSurface();
		~MeshExtractOuterSurface();

		// accessors
		/**
		* Sets the volume mesh
		*/
		void SetVolumeMesh (vtkPolyData* volumeMesh) {_volumeMesh = volumeMesh;}
		/**
		* Gets the surface mesh
		*/
		vtkPolyData* GetSurfaceMesh () {return _surfaceMesh;}

		// member functions
		/**
		* Extracts the outer surface from a volume mesh
		*/
		void Update();
		
	
	protected:
	
		// data members from accessors
		vtkPolyData* _volumeMesh;
		vtkPolyData* _surfaceMesh;		
		

	};
}


// Siegrist's mesh type
class FEMAssignerMesh{

public:
		
	//A surface triangle of three node Ids
	struct Triangle {
		int a;
		int b;
		int c;
	};

	// An element with its node Ids
	struct Element {
		std::vector<int> triangleNodes;
		//std::vector<int> intermediateNodes;
	};

	//The surface triangles
	std::vector<FEMAssignerMesh::Triangle> triangleVector;
	
	//A map that maps the element Ids to the mesh::Mesh::Element data
	std::map<int, Element> elements;
};


struct ElementTriangle {
	int elementId;
	FEMAssignerMesh::Triangle triangle;
};

#endif //MESHEXTRACTOUTERSURFACE_H