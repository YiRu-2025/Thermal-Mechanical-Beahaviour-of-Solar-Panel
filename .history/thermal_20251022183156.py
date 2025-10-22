# Thermal analysis module 
# heat transfer and temperature distribution calculation by finite element method

import modelget
import gmsh
import numpy as np  

# step 1: element conductivity matrix and load vector(heat flux) calculation

# Thermal conductivity:
al_tc = 200.0 # W/m·K
si_tc = 150.0 # W/m·K
glass_tc = 1.0 # W/m·K

# Initialize Gmsh and load the mesh

gmsh.open("solar_panel.msh")

# Get all nodes
#   node_tags: node numeration
#   node_coords: node cordinates [x,y]
#   elements_al:  elements (defined by nodes) in the aluminum frame
#   elements_si:  elements (defined by nodes) in the cells
#   elements_glass: elements (defined by nodes) in the glass cover

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_coords = np.array(node_coords).reshape(-1, 3)
print("Nodes = ", node_coords.shape)

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
elements_surface_al = get_elements_of_physical_group(2, 1)
elements_surface_si = get_elements_of_physical_group(2, 2)
elements_surface_gl = get_elements_of_physical_group(2, 3)

print(elements_surface_al.shape, elements_surface_si.shape, elements_surface_gl.shape)

# define cst function to get s and p
def jacobian(coord):
  # compute the jacobian value for a triangle
  # coord is a (2x3) array
  # [[x1 x2 x3]
  #   y1 y2 y3]]

  jacobian = (coord[0,0]-coord[0,2]) * (coord[1,1]-coord[1,2]) - (coord[0,1] - coord[0,2]) * (coord[1,0] - coord[1,2])
  # it's the determinant of the jacobian
  return abs(jacobian)

def b_matrix(coord):
    B = np.zeros((2, 3))
    B[0, :] = np.roll(coord[1, :], -1) - coord[1, :]
    B[1, :] = coord[0, :] - np.roll(coord[0, :], -1)
    B /= jacobian(coord)
    return B

def CST(coord, conductivity, flux):
  B = b_matrix(coord)
  # multiply by the conductivity value and by the area of the element and also by the thickness!
  # Compute transpose(B).B using numpy.dot
  s = np.dot(B.T, B) * conductivity * jacobian(coord)

  p =np.zeros((3))

  if coord[1,0] == 0.06 and coord[1,1] == 0.06:
    #print('Nodes 1 and 2 are on the top edge')
    distance = np.absolute(coord[0,0] - coord[0,1])
    p[0] = flux*distance/2
    p[1] = flux*distance/2
  elif coord[1,0] == 0.06 and coord[1,2] == 0.06:
    #print('Nodes 1 and 3 are on the top edge')
    distance = np.absolute(coord[0,0] - coord[0,2])
    p[0] = flux*distance/2
    p[2] = flux*distance/2
  elif coord[1,1] == 0.06 and coord[1,2] == 0.06:
    #print('Nodes 2 and 3 are on the top edge')
    distance = np.absolute(coord[0,0] - coord[0,1])
    p[1] = flux*distance/2
    p[2] = flux*distance/2

  return s,p

# step 2 Assembly global stiffness matrix and load vector
# Here we start with the FE process
# In hands we have:
# the nodes coordinates: node_coords
# the connectivity for each material: elements_surface_1, elements_surface_2, elements_surface_3

# Remove all third coordinates (should be zero's!)
nodes = np.delete(node_coords, 2, 1)
nn = nodes.shape[0]
print("Nodes = ", nodes.shape)
ne_1 = elements_surface_al.shape[0]
ne_2 = elements_surface_si.shape[0]
ne_3 = elements_surface_gl.shape[0]
#print("Aluminium support = ", elements_surface_1.shape)
#print("Cells = ", elements_surface_2.shape)
#print("Glass cover = ", elements_surface_3.shape)


# Modify the calculate_element_stiffness function as provided earlier

def assemble(element, s, p):
    e_index = np.array(element)-1   # start from 0

    for i in range(len(element)):
        # Assemble the element stiffness matrix into the global stiffness matrix
        for j in range(len(element)):
            K[e_index[i], e_index[j]] += s[i, j]
        #if np.count_nonzero(p) > 0:
        r[e_index[i]] += p[i]
    return K, r


K = np.zeros((nn, nn))
r = np.zeros(nn)

# Loop over all the elements and compute the elementary quantities using the "CST" function
# Accumulate contributions into the global stiffness matrix K and the global right-hand side vector r
for element in elements_surface_al:
    coord = np.transpose(np.vstack((nodes[int(element[0]-1)], nodes[int(element[1]-1)], nodes[int(element[2]-1)])))
    s, p = CST(coord=coord, conductivity=al_tc, flux=1000.0)
    K, r = assemble(element, s, p)

for element in elements_surface_si:
    coord = np.transpose(np.vstack((nodes[int(element[0]-1)], nodes[int(element[1]-1)], nodes[int(element[2]-1)])))
    s, p = CST(coord=coord, conductivity=si_tc, flux=1000.0)
    K, r = assemble(element, s, p)

for element in elements_surface_gl:
    coord = np.transpose(np.vstack((nodes[int(element[0]-1)], nodes[int(element[1]-1)], nodes[int(element[2]-1)])))
    s, p = CST(coord=coord, conductivity=glass_tc, flux=1000.0)
    K, r = assemble(element, s, p)

print("Global Stiffness Matrix shape:", K.shape)
#print(K)

print("Global Right-hand Side Vector shape:", r.shape)
#print(r)

# step 3 apply boundary conditions and solve the linear system
# method 1: delete the rows and columns corresponding to the fixed nodes in K and r
# The equation need to solve is K * T = r
# We already get K and r, now we need to solve the temperature T

# The thermal BC:
# top edge : heat flux = 1000 W/(m^2) --Neumann BC
# left edge :
# right edge :
# bottom edge : T = 0 --Dirichlet BC

### Imposing Dirichlet BC by Deleting the rows and columns of the K and r to simplify the matrix
# 1st: need to find all the bottom nodes -- bottom_nodes
# 2nd: delete them in K and r at the same time!! Get the dense K and r from the corse one
# 3rd: compute T
# 4th: add the delete nodes back in T with value 0

# step 1: find all the bottom nodes
bottom_nodes = np.where(nodes[:, 1] == 0)[0]  # Use [0] to get the array of indices
print("Indices of value in array:", bottom_nodes)

# step 2: delete all the bottom nodes in K and r meanwhile
K_deleto = np.delete(np.delete(K, bottom_nodes, axis=0), bottom_nodes, axis=1)  # Delete corresponding rows of K

r_deleto = np.delete(r, bottom_nodes, axis=0)  # Delete corresponding rows of r

print("K_deleto size: ",K_deleto.shape)
# checking if K and r deleted well
#print(K.shape)
#print(r.shape)

# step 3: compute T
T = np.linalg.solve(K_deleto, r_deleto)
print(T.shape)

# step 4: Add the deleted rows to the final T vector
for i in bottom_nodes:
 T = np.insert(T, i, 0, axis=0)
#print(T)

# recover the K_deleto and r_deleto to avoid running times problem
K_deleto = K
r_deleto = r
print("K_recoery size: ",K_deleto.shape)

# method 2: replace 0 at the corresponding at the corresponding positions and 1 at diagonal positions
# This code works well but not good for large matrix solving efficiency
'''
import matplotlib.pyplot as plt
import matplotlib.tri as tri
j=0
aa=np.zeros(40)
####Implementation of the Dirichlet BC: @ y=0 : T(y=0)=T0=0
for element in elements_surface_al:
  coord = np.transpose(np.vstack((nodes[int(element[0]-1)], nodes[int(element[1]-1)], nodes[int(element[2]-1)])))
  for i in [0,1,2]:
   if coord[1,i] == 0 :   # finidng the nodes on the y=0 line
    aa[j]=element[i]-1
    j=j+1
    K[int(element[i]-1),:]=0
    K[:,int(element[i]-1)]=0
    r[int(element[i]-1)]=0
    K[int(element[i]-1),int(element[i]-1)]=1

####Solution
T = np.linalg.solve(K, r)

'''

# structural flux calculation
# Here I calculate the flux of each element by using the centoriod since we consider the value as a constant
# in order to see how the heat flow moves inside the panel

def calculate_flux_and_centroid(elements, nodes, T, material_tc):
    elem_num = len(elements)
    flux = np.zeros((2, elem_num))
    centroid = np.zeros((2, elem_num))

    for i, element in enumerate(elements):
        coord = np.transpose(np.vstack((nodes[int(element[0]-1)], nodes[int(element[1]-1)], nodes[int(element[2]-1)])))
        T_elem = [T[int(element[0]-1)], T[int(element[1]-1)], T[int(element[2]-1)]]
        B = b_matrix(coord)

        flux[:, i] = np.dot(B, np.transpose(T_elem)) * (-material_tc)
        centroid[0, i] = np.mean(coord[0, :])
        centroid[1, i] = np.mean(coord[1, :])

    return flux, centroid

flux1, centroid1 = calculate_flux_and_centroid(elements_surface_al, nodes, T, al_tc)
flux2, centroid2 = calculate_flux_and_centroid(elements_surface_si, nodes, T, si_tc)
flux3, centroid3 = calculate_flux_and_centroid(elements_surface_gl, nodes, T, glass_tc)

# Combine flux and centroid arrays
flux = np.hstack((flux1, flux2, flux3))
centroid = np.hstack((centroid1, centroid2, centroid3))

# step 4: post-process and visualize the results for thermal analysis

# temperature distribution plotting
# using matplotlib's tricontour and tricontourf functions
# need to triangutalte the coordinates to plot the temperature results
import matplotlib.pyplot as plt
import matplotlib.tri as tri

triang = tri.Triangulation(node_coords[:, 0], node_coords[:, 1])

plt.figure(figsize=(18,8))
plt.tricontour(triang, T, levels=20, colors='g', linewidths=0.5) # Plot the equivalent temperature line
plt.tricontourf(triang, T, levels=20, cmap='inferno') # Plot the temperature area

plt.gca().set_aspect('equal', adjustable='box')

plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.title('Temperature')
plt.colorbar(label='Temperature (°C)')

plt.show()

# heat flow plotting
plt.quiver(centroid[0,:], centroid[1,:], flux[0,:], flux[1,:])
plt.show()