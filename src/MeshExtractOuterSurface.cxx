/*
 * ISTB - University of Bern, Serena Bonaretti
 */

#include <MeshExtractOuterSurface.h>

#include <PointWriterXyz.h>

#include <map>
#include <vector>
#include <set>

#include <vtkCellArray.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkSTLWriter.h>

using namespace std;
using namespace points;


// it sorts the nodes of each triangle of the tetra in an ascending order
void insertTriangle(ElementTriangle* etriangle, map<int,vector<ElementTriangle> >* triangleMap, int a, int b, int c) {
	
	if (a > b)
		swap(a, b);
	if (b > c)
		swap(b, c);
	if (b < a)
		swap(a, b);

	etriangle->triangle.a = a;
	etriangle->triangle.b = b;
	etriangle->triangle.c = c;

	vector<ElementTriangle>* etriangles = &((*triangleMap)[a]);
	etriangles->push_back(*etriangle);
}

// it compares if the sorted nodes of two triangles are the same
bool same(FEMAssignerMesh::Triangle* t1, FEMAssignerMesh::Triangle* t2) {
	return (t1->a == t2->a && t1->b == t2->b && t1->c == t2->c);
}



namespace mesh {

	// constructor
	MeshExtractOuterSurface::MeshExtractOuterSurface(){

		_surfaceMesh = vtkPolyData::New();

	}


	// destructor
	MeshExtractOuterSurface::~MeshExtractOuterSurface(){
	
		_surfaceMesh->Delete();

	}


	// member function
	void MeshExtractOuterSurface::Update(){

		FEMAssignerMesh* mesh = new FEMAssignerMesh;
		FEMAssignerMesh::Element element;
		vtkCellArray* polys = vtkCellArray::New();
		vtkIdList* idList;
		vtkGenericCell* cell = vtkGenericCell::New();
		
		// insert first four nodes in Element of FEMAssignerMesh from vtkPolyData
		for (int i=0; i<_volumeMesh->GetNumberOfCells(); i++){

			// get the cell
			_volumeMesh->GetCell(i, cell);
			
			// get the 4 tetra vertices
			idList = cell->GetPointIds();
			
			// set them in mesh
			for (int y=0; y<4; y++) {
				element.triangleNodes.push_back(idList->GetId(y));
			}
			mesh->elements.insert(make_pair(i, element));
			element.triangleNodes.clear();
		}

		// maps that contains all sorted triangle whose first node is the same
		std::map<int, vector<ElementTriangle> > allTriangles; // map with first triangle node as index

		for (map<int, FEMAssignerMesh::Element>::iterator it = mesh->elements.begin(); it
			!= mesh->elements.end(); ++it) {
		int elementId = (*it).first;

		std::vector<int>* tetrahedron = &(*it).second.triangleNodes;

		// inserts all triangle permutations
		ElementTriangle etriangle;
		etriangle.elementId = elementId;
		
		insertTriangle(&etriangle, &allTriangles, tetrahedron->at(0),
				tetrahedron->at(1), tetrahedron->at(2)); // 012
		insertTriangle(&etriangle, &allTriangles, tetrahedron->at(0),
				tetrahedron->at(1), tetrahedron->at(3)); // 013
		insertTriangle(&etriangle, &allTriangles, tetrahedron->at(0),
				tetrahedron->at(2), tetrahedron->at(3)); // 023
		insertTriangle(&etriangle, &allTriangles, tetrahedron->at(1),
				tetrahedron->at(2), tetrahedron->at(3)); // 123
		}

		// looks for non-shared faces
		std::set <int> nodeIds;
		for (std::map<int, std::vector<ElementTriangle> >::iterator it = allTriangles.begin(); it
			!= allTriangles.end(); ++it) {
				std::vector<ElementTriangle>* etriangles = &(*it).second;
				for (int i = 0; i < etriangles->size(); i++) {
					if (etriangles->at(i).elementId == -1)
						continue;
					bool unique = true;
					for (int j = i + 1; j < etriangles->size(); j++) {
						if (same(&etriangles->at(i).triangle,
								&etriangles->at(j).triangle)) {
							unique = false;
							etriangles->at(j).elementId = -1;
						}
					}
					if (unique) {
						
						// faces back to vtkPolyData for the surface mesh
						polys->InsertNextCell(3);
						polys->InsertCellPoint(etriangles->at(i).triangle.a);
						polys->InsertCellPoint(etriangles->at(i).triangle.b);
						polys->InsertCellPoint(etriangles->at(i).triangle.c);
												
						// surface nodes without replication
						std::set<int>::iterator it; 
						if (nodeIds.count(etriangles->at(i).triangle.a) == 0){
							nodeIds.insert(etriangles->at(i).triangle.a);
						}
						if (nodeIds.count(etriangles->at(i).triangle.b) == 0){
							nodeIds.insert(etriangles->at(i).triangle.b);
						}
						if (nodeIds.count(etriangles->at(i).triangle.c) == 0){
							nodeIds.insert(etriangles->at(i).triangle.c);
						}
					}
				}
		}

		// sets faces
		_surfaceMesh->SetPolys(polys);

		// sets nodes
		vtkPoints* points = vtkPoints::New();
		std::set<int>::iterator it;
		for (it = nodeIds.begin(); it != nodeIds.end(); ++it){
			double pt[3];
			_volumeMesh->GetPoint(*it, pt);
			points->InsertNextPoint(pt);
		}
		_surfaceMesh ->SetPoints(points);

		//std::cout << "number of total points: " << _volumeMesh->GetNumberOfPoints() << std::endl; 
		//std::cout << "number of surface points: " << _surfaceMesh->GetNumberOfPoints() << std::endl; 
		

		// cleaning up
		cell->Delete();
		delete mesh;
		polys->Delete();
		points->Delete();
				
	}
}
	