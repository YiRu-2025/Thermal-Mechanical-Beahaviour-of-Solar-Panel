# The same procedures as thermal analysis. 
# Need to modify the CST function to mechanical case and expand the K and r to twice of the node size.
import numpy as np
import thermal
#step 1 elementary stiffness Ke and elementary force vector fe

# Here we define the CST_mechanic function to obtain Ke and fe
#
# input: nodes_coordinate
#     materail parameter list
#     nodal temperature vector T
#
# output: Ke(6x6), fe(6x1)
#
# inside the function:
# Ke = det(J) * thickness * B.transpose * C * B   , independent to T
# fe = det(J) * thickness * B.transpose * stress
#   = det(J) * thickness * B.transpose * (C * (-alpha * avg_T * I)) ,
# where B becomes a 3x6 matrix need to rewrite,
#     alpha is the thermal expansion coefficient,
#     avg_T is the average temperature of 3 nodes
#     C is the 3x3 modulus matrix
def compute_C(E,v):
  C = E/(1.0+v)/(1.0-2.0*v) * np.array([[ 1.0-v,    v,   0.0],
								        [   v,  1.0-v,   0.0],
								        [  0.0,   0.0,   0.5-v]])
  return C

def material_properties(material_type):
    properties = {}

    if material_type == "Aluminum":
        properties['E'] = 68.9e9  # Young's modulus (Pa)
        properties['v'] = 0.33  # Poisson's ratio
        properties['a'] = 23.1e-6  # Thermal expansion coefficient (K^(-1))
        properties['rho'] = 2700  # Density (kg/m^3)
        properties['yield_stress'] = 276e6  # Elastic limit (Pa)
        properties['C'] = compute_C(properties['E'],properties['v'])
    elif material_type == "Silicon":
        properties['E'] = 74.8e9
        properties['v'] = 0.19
        properties['a'] = 0.55e-6
        properties['rho'] = 2650
        properties['yield_stress'] = 100e6
        properties['C'] = compute_C(properties['E'],properties['v'])
    elif material_type == "Glass":
        properties['E'] = 72.0e9
        properties['v'] = 0.22
        properties['a'] = 8.0e-6
        properties['rho'] = 2200
        properties['yield_stress'] = 30e6
        properties['C'] = compute_C(properties['E'],properties['v'])

    else:
        raise ValueError("Invalid material type. Supported types: 'Aluminum', 'Silicon', 'Glass'.")

    return properties


# Example usage:
aluminum_properties = material_properties("Aluminum")
silicon_properties = material_properties("Silicon")
glass_properties = material_properties("Glass")


def B_matrix(coord):
  b = np.zeros((3, 6))
  x1, y1 = coord[0, 0], coord[1, 0]
  x2, y2 = coord[0, 1], coord[1, 1]
  x3, y3 = coord[0, 2], coord[1, 2]

  dN1_dx = (y2 - y3)
  dN2_dx = (y3 - y1)
  dN3_dx = (y1 - y2)

  dN1_dy = (x3 - x2)
  dN2_dy = (x1 - x3)
  dN3_dy = (x2 - x1)

  # B matrix
  b = np.array([
        [dN1_dx,  0,    dN2_dx,  0,   dN3_dx,  0],
        [0,    dN1_dy,  0,    dN2_dy, 0,    dN3_dy],
        [dN1_dy, dN1_dx,  dN2_dy,  dN2_dx, dN3_dy, dN3_dx]
      ])
  b /= thermal.jacobian(coord)
  return b

def compute_thermal_strain(element,a):  # Thermal strain
  elem_index = np.array(element) - 1

  avg_T = np.sum(thermal.T[elem_index])/3
  strain_t = a * (avg_T) * np.array([1,1,0])
  return strain_t

def CST_mechanic(coord,C):

  b = B_matrix(coord)

  Ke = (thermal.jacobian(coord)*0.5) * np.dot(np.dot(np.transpose(b),C),b) # 6x6 matrix

  fe = np.zeros((6))
  fe = (thermal.jacobian(coord)*0.5) * np.dot(np.transpose(b),C)

  return Ke, fe

##Step 2: Assemble global stiffness and force vector
## Here I assemble the Ke and re to global stiffness and force vector combining step 1
# only need to modify the size and index of K and r based on previous thermal case

def assemble_mechanic(element, ke, fe):
    e_index = np.array(element)-1   # start from 0

    for i in range(len(element)):
        # Assemble the element stiffness matrix into the global stiffness matrix
        for j in range(len(element)):
            K[2 * e_index + i % 2, 2 * e_index + j % 2] += ke[i, j]

        f[2 * e_index + i % 2] += fe[i]

    return K, f


# Initialize global stiffness matrix and right-hand side vector
nn = thermal.nn  # number of nodes
K = np.zeros((nn*2, nn*2))
f = np.zeros(nn*2)

# Loop over all the elements and compute the elementary quantities using the "CST" function
# Accumulate contributions into the global stiffness matrix K and the global right-hand side vector r
nodes = thermal.nodes
for element in thermal.elements_surface_al:
    coord = np.transpose(np.vstack((nodes[int(element[0] - 1)], nodes[int(element[1] - 1)], nodes[int(element[2] - 1)])))
    C = aluminum_properties['C']
    ke,re = CST_mechanic(coord,C)
    K, f = assemble_mechanic(element, ke, re)

# Repeat the process for elements_surface_2 and elements_surface_3
for element in thermal.elements_surface_si:
    coord = np.transpose(np.vstack((nodes[int(element[0] - 1)], nodes[int(element[1] - 1)], nodes[int(element[2] - 1)])))
    C = silicon_properties['C']
    ke,re = CST_mechanic(coord,C)
    K, f = assemble_mechanic(element, ke, re)

for element in thermal.elements_surface_gl:
    coord = np.transpose(np.vstack((nodes[int(element[0] - 1)], nodes[int(element[1] - 1)], nodes[int(element[2] - 1)])))
    C = glass_properties['C']
    ke,re = CST_mechanic(coord,C)
    K, f = assemble_mechanic(element, ke, re)

# Now we have the global stiffness matrix K and the global force vector r

# step 3: Apply boundary conditions and solve the system
## Imposing Dirichlet BC: Second method- Deleting the associated rows and columns of the K and the associated rows of the f vector

# 1- ### Implementation of Dirichlet BC: @ y=0 : u(y=0)=v(y=0)=0
# finding the corresponding nodes to apply Dirichlet BC at the bottom line on y=0
#nodes_bottom = np.zeros(11)    ## nodes at the bottom of the mesh (at y=0 )

nodes_bottom = np.where(nodes[:,1] == 0)[0]        ## finding the nodes at the bottom of the mesh : @ y=0
nodes_left = np.where(nodes[:,0] == 0)[0]          ## finding the nodes at the left side of the mesh: @ x=0
nodes_center = np.where(nodes[:,0] == 0.1)[0]         ## finding the nodes on the center line of the mesh: @ x=0.1

# Exclude the shared nodes at (0,0) and (0.1,0) from nodes_left and nodes_center
nodes_left = np.setdiff1d(nodes_left, nodes_bottom)
nodes_center = np.setdiff1d(nodes_center, nodes_bottom)

print("Node numbers of the nodes at the bottom line:", nodes_bottom)
print("Node numbers of the nodes at the left side:", nodes_left)
print("Node numbers of the nodes on the centeral line:", nodes_center)


rows_to_remove = np.zeros(2*(len(nodes_bottom))+ len(nodes_left)+ len(nodes_center), dtype=int)
cols_to_remove = rows_to_remove

# 2- Removing the Corresponding rows and columns of the K and r
i=0
for i, node_b in enumerate(nodes_bottom):
  rows_to_remove[2*i] = 2*int(node_b)
  rows_to_remove[2*i+1] = 2*int(node_b)+1
  cols_to_remove[2*i] = 2*int(node_b)
  cols_to_remove[2*i+1] = 2*int(node_b)+1
for i, node_l in enumerate(nodes_left):   # here we need to just remove the components of the rows associated to the u and leave the y ones as the BC here is only u(@x=0)=0 and v is free to move
  rows_to_remove[2*len(nodes_bottom)+i] = 2*int(node_l)
  cols_to_remove[2*len(nodes_bottom)+i] = 2*int(node_l)

for i, node_c in enumerate(nodes_center):  # here we need to just remove the components of the rows associated to the u and leave the y ones as the BC here is only u(@x=0.1)=0 and v is free to move
  rows_to_remove[2*len(nodes_bottom)+len(nodes_left)+i] = 2*int(node_c)
  cols_to_remove[2*len(nodes_bottom)+len(nodes_left)+i] = 2*int(node_c)


# Remove rows and columns from K
K_deleto = np.delete(np.delete(K, rows_to_remove, axis=0), cols_to_remove, axis=1)
# Remove rows from r
f_deleto = np.delete(f, rows_to_remove, axis=0)

print("K_deleto size: ",K_deleto.shape)
#### Solution

# 1- Solving the reduced system of equations
U = np.linalg.solve(K_deleto, f_deleto)
print("U size without 0:",U.shape)
# 2- Adding the deleted rows to the final T vector in the associated positions of the nodes
for i in rows_to_remove:
 U = np.insert(U, i, 0, axis=0)
print("U total size:",U.shape)

# recover the K_deleto and r_deleto to avoid running times problem
K_deleto = K
f_deleto = f
print("K_recoery size: ",K_deleto.shape)

# method 2: replacing
'''
def apply_fixed_bc(K, f, idx):
  K[idx,:]=0
  K[:,idx]=0
  K[idx+1,:]=0
  K[:,idx+1]=0

  f[idx]=0
  f[idx+1]=0

  K[idx,idx]=1
  K[idx+1,idx+1]=1

def apply_slide_bc(K, f, idx):
  K[idx,:]=0
  K[:,idx]=0

  f[idx]=0
  K[idx,idx]=1

idx = 0
####Implementation of Dirichlet BC: @ y=0 : u(y=0)=v(y=0)=0
for element in elements_surface_al:
  coord = np.transpose(np.vstack((nodes[int(element[0]-1)], nodes[int(element[1]-1)], nodes[int(element[2]-1)])))
  for i in [0,1,2]:
    idx = 2*int(element[i]-1)
    if coord[1,i] == 0 :   # bottom
      apply_fixed_bc(K, f, idx)
    elif coord[0,i] == 0 : # left
      apply_slide_bc(K, f, idx)
    elif coord[0,i] == 0 :   # right
      apply_slide_bc(K, f, idx)


####Implementation of Dirichlet BC: @ x=0 : u(x=0)=0
for element in elements_surface_si:
  coord = np.transpose(np.vstack((nodes[int(element[0]-1)], nodes[int(element[1]-1)], nodes[int(element[2]-1)])))
  for i in [0,1,2]:
    if coord[0,i] == 0 :   # right
      idx = 2*int(element[i]-1)
      apply_slide_bc(K, f, idx)


####Implementation of Dirichlet BC: @ x=0.1 : u(x=0.1)=0
for element in elements_surface_gl:
  coord = np.transpose(np.vstack((nodes[int(element[0]-1)], nodes[int(element[1]-1)], nodes[int(element[2]-1)])))
  for i in [0,1,2]:
   if coord[0,i] == 0.1 :   # right
    idx = 2*int(element[i]-1)
    apply_slide_bc(K, f, idx)
# Solution
d = np.linalg.solve(K, f)
num_zeros = np.count_nonzero(d == 0)
print("Number of zeros in d:", num_zeros)

'''

# step 4: thermal-mechancial behavior analysis: stress, strain, displacement
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# Extracting u and v components from U
u = U[::2]  # Every other element starting from index 0
v = U[1::2]  # Every other element starting from index 1
print('max u=', max(u))
print('max v=', max(v))
# Plotting the structure with displacement vectors
plt.figure(figsize=(10, 6))

# Plot the original structure (nodes without displacement)
node_coords = thermal.nodes_coords
plt.plot(node_coords[:, 0], node_coords[:, 1], 'o', color='blue', markersize=5, label='Original Structure')

# Plot the deformed structure using displacement vectors
plt.quiver(node_coords[:, 0], node_coords[:, 1], u, v, scale=1, scale_units='xy', angles='xy', color='red', label='Displacement Vectors')

plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.title('Deformed Structure with Displacement Vectors')
plt.legend()
plt.grid(True)

# Compute the magnitude of the displacement at each node
displacement_magnitude = np.sqrt(u**2 + v**2)

# Triangulation for contour plot
triang = tri.Triangulation(node_coords[:, 0], node_coords[:, 1])

# Plotting the vector plot for the total displacement U
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_aspect('equal')

# Using quiver to plot arrows indicating the displacement direction
q = ax.quiver(node_coords[:, 0], node_coords[:, 1], u, v, scale=50, scale_units='xy', angles='xy', color='blue', label='Displacement')

# Plotting the contour plot for the displacement magnitude
tcf = ax.tricontourf(triang, displacement_magnitude, levels=20, cmap='viridis', alpha=0.7, label='Displacement Magnitude')

# Adding a colorbar for the displacement magnitude
fig.colorbar(tcf, label='Displacement Magnitude (m)')

ax.set_title('Vector Plot of Total Displacement with Displacement Magnitude Contour')
ax.set_xlabel('X coordinate')
ax.set_ylabel('Y coordinate')

plt.legend()
plt.show()

####Solution
# Extracting u and v components from U
# u = U[::2]  # Every other element starting from index 0
# v = U[1::2]  # Every other element starting from index 1


####Contour Plot
# plot u
fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tcf = ax1.tricontourf(triang,u)
fig1.colorbar(tcf)
ax1.set_title('Contour Plot of Displacement in X direction')

# plot v
fig2, ax2 = plt.subplots()
ax2.set_aspect('equal')
tcf = ax2.tricontourf(triang, v)
fig2.colorbar(tcf)
ax2.set_title('Contour Plot of Displacement in Y direction')

# Strain: ε = [B]·[U] # matrix order : 3x1 = 3x6 · 6x1
# Stress: σ = C : (ε - εt)
'''
Put my main idea here
### step 1: get the elementary total strain and elementary stress
for element in each surface

  thermal_expansion = surface_properties['a']
  C = surface_properties['C']

  ## compute elementary strain and stress -- can be written as a function

### step 2: assemble strain_elem amd stress_elem to global one

### step 3: get the von_mises equivalent stress
'''
node_tags = thermal.node_tags
u = np.transpose(np.array([U[:len(node_tags)],U[len(node_tags):]]))
strain_elem = np.zeros((3, 1))
stress_elem = np.zeros((3, 3))
def compute_element_stress_strain(element,C,a):

  e_index = np.array(element) - 1 # index starts from 0
  coord = np.transpose(np.vstack((nodes[e_index[0]], nodes[e_index[1]], nodes[e_index[2]])))

  # get the element displacement vector 6x1
  U_elem = np.array([u[e_index[0], 0], u[e_index[0], 1], u[e_index[1], 0], u[e_index[1], 1], u[e_index[2], 0], u[e_index[2], 1]])

  strain_elem = np.dot(B_matrix(coord),U_elem) # 3x1 vector
  # strain_elem[2] = 2 * strain_elem[2] # 2 * (εxy)
  # It's unnecessary when converting the strain tensor to its vector form
  # since I don't keep the half factor in the derivative of the shape function: B
  ## compute thermal strain
  strain_thermal = compute_thermal_strain(element, a)
  stress_elem = np.dot(C, strain_elem - strain_thermal)

  return strain_elem, stress_elem

## assemble strain&stress
def assemble_strain_stress(element,C,a):
  e_index = np.array(element) - 1
  strain_elem, stress_elem = compute_element_stress_strain(element,C,a)
  for node in e_index:
    for i in range(3):
      strain[node][i] += strain_elem[i]
      stress[node][i] += stress_elem[i]

  return stress, strain

stress = np.zeros((len(node_tags),3))
strain = np.zeros((len(node_tags),3))

## loop and assemble
for element in thermal.elements_surface_al:
  a = aluminum_properties['a']
  C = aluminum_properties['C']
  assemble_strain_stress(element,C,a)
for element in thermal.elements_surface_si:
  a = silicon_properties['a']
  C = silicon_properties['C']
  assemble_strain_stress(element,C,a)
for element in thermal.elements_surface_gl:
  a = glass_properties['a']
  C = glass_properties['C']
  assemble_strain_stress(element,C,a)
# Now we get total stain and stress

import matplotlib.pyplot as plt
import matplotlib.tri as tri

plt.figure(figsize=(18,8))
plt.tricontour(triang, abs(stress[:,0]), levels=20, colors='k', linewidths=0.5)
plt.tricontourf(triang, abs(stress[:,0]), levels=20, cmap='jet')

plt.plot(node_coords[:, 0], node_coords[:, 1], 'o', color='red', markersize=1)

plt.gca().set_aspect('equal', adjustable='box')

plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.title('σxx')
plt.colorbar(label='Stress (Pa)')


plt.figure(figsize=(18,8))
plt.tricontour(triang, abs(stress[:,1]), levels=20, colors='k', linewidths=0.5)
plt.tricontourf(triang, abs(stress[:,1]), levels=20, cmap='jet')

plt.plot(node_coords[:, 0], node_coords[:, 1], 'o', color='red', markersize=1)

plt.gca().set_aspect('equal', adjustable='box')

plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.title('σyy')
plt.colorbar(label='Stress (Pa)')


plt.show()