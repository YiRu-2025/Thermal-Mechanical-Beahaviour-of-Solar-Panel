import gmsh
import numpy as np
import matplotlib.pyplot as plt


# Initialize Gmsh
gmsh.initialize()

# Create a new model
gmsh.model.add("model")

######## GEOMETRY DEFINITION

# 7-----6---------------5
# |     |   material A  |
# |     8---------------4
# |     |   material B  |            
# |     9---------------3
# |  material C         |                   
# 1---------------------2

lc1 = 0.005
lc2 = 0.01  # characteristic lengths

# Define points
gmsh.model.geo.addPoint(0, 0, 0, lc2, 1)
gmsh.model.geo.addPoint(0.1, 0, 0, lc2, 2)
gmsh.model.geo.addPoint(0.1, 0.04, 0, lc1, 3)
gmsh.model.geo.addPoint(0.1, 0.05, 0, lc1, 4)
gmsh.model.geo.addPoint(0.1, 0.06, 0, lc1, 5)
gmsh.model.geo.addPoint(0.04, 0.06, 0, lc1, 6)
gmsh.model.geo.addPoint(0, 0.06, 0, lc2, 7)
gmsh.model.geo.addPoint(0.04, 0.05, 0, lc1, 8)
gmsh.model.geo.addPoint(0.04, 0.04, 0, lc1, 9)


# Define lines
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 5, 4)
gmsh.model.geo.addLine(5, 6, 5)
gmsh.model.geo.addLine(6, 7, 6)
gmsh.model.geo.addLine(7, 1, 7)
gmsh.model.geo.addLine(6, 8, 8)
gmsh.model.geo.addLine(8, 4, 9)
gmsh.model.geo.addLine(8, 9, 10)
gmsh.model.geo.addLine(9, 3, 11)


# Define the loop and surface for aluminium panel
gmsh.model.geo.addCurveLoop([1, 2, -11, -10, -8, 6, 7], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

# Define the loop and surface for cells
gmsh.model.geo.addCurveLoop([3, -9, 10, 11], 2)
gmsh.model.geo.addPlaneSurface([2], 2)

# Define the loop and surface for cover
gmsh.model.geo.addCurveLoop([8, 9, 4, 5], 3)
gmsh.model.geo.addPlaneSurface([3], 3)

# Synchronize the CAD kernel with the Gmsh model
gmsh.model.geo.synchronize()

# Define two physical surfaces
gmsh.model.addPhysicalGroup(2, [1], 1)
gmsh.model.setPhysicalName(2, 1, "Physical Surface 1")

gmsh.model.addPhysicalGroup(2, [2], 2)
gmsh.model.setPhysicalName(2, 2, "Physical Surface 2")

gmsh.model.addPhysicalGroup(2, [3], 3)
gmsh.model.setPhysicalName(2, 3, "Physical Surface 3")

# Generate the mesh
gmsh.model.mesh.generate(2)

# Save and close
gmsh.write("solar_panel.msh")

# Initialize Gmsh and load the mesh
gmsh.initialize()
gmsh.open("solar_panel.msh")

# Get all nodes
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_coords = np.array(node_coords).reshape(-1, 3)

# Function to get elements of a physical group
def get_elements_of_physical_group(dim, tag):
    entity_tags = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    all_elements = []
    for entity_tag in entity_tags:
        element_types, element_tags, node_tags = gmsh.model.mesh.getElements(dim, entity_tag)
        for i in range(len(element_types)):
            all_elements.append(np.array(node_tags[i]).reshape(-1, 3))
    return np.concatenate(all_elements) if all_elements else np.array([])

# Get elements for each physical surface
elements_surface_1 = get_elements_of_physical_group(2, 1)
elements_surface_2 = get_elements_of_physical_group(2, 2)
elements_surface_3 = get_elements_of_physical_group(2, 3)

# Plot the mesh for each physical surface separately
plt.figure(figsize=(8, 8))

# Plot for Physical Surface 1
if elements_surface_1.size > 0:
    plt.triplot(node_coords[:, 0], node_coords[:, 1], elements_surface_1 - 1, color='blue', label='Aluminium support')

# Plot for Physical Surface 2
if elements_surface_2.size > 0:
    plt.triplot(node_coords[:, 0], node_coords[:, 1], elements_surface_2 - 1, color='red', label='Cells')

# Plot for Physical Surface 3
if elements_surface_3.size > 0:
    plt.triplot(node_coords[:, 0], node_coords[:, 1], elements_surface_3 - 1, color='green', label='Glass cover')

plt.gca().set_aspect('equal')
plt.title('Solar panel mesh with three physical surfaces')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()

# Finalize Gmsh
gmsh.finalize()


