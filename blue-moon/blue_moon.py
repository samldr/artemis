import gmsh
import math

gmsh.initialize()

######## MODEL PARAMETERS ########

height = 16 # overall lander height
radius = 3 # fuselage radius
tank_radius = 1.25 # radial tank radius
leg_radius = 0.2 # landing leg radius
boundary_radius = 35 # boundary cylinder radius
boundary_height = 15 # bboundary cylinder height
tolerance = 0.05 # empty space between surfaces of different physical groups

######## MESH PARAMETERS ########

meshsize_lowerfuselage = 0.2 * radius # upper fuselage
meshsize_upperfuselage = 0.2 * (radius - 0.5) # lower fuselage
meshsize_tanks = 0.2 * tank_radius # radially mounted tanks
meshsize_landinglegs = 0.3 * leg_radius # landing legs
meshsize_space = 0.5 # top and side boundaries (space)
meshsize_ground = 0.5 # lunar surface boundary

# function to set the mesh size of surfaces from its physical group tag
def setMeshSize(physical_group, mesh_size):
    
    gmsh.model.occ.synchronize()

    # gets the lines that make up each physical surface, and then the points that make up each of those lines
    points = []
    for tag in gmsh.model.getEntitiesForPhysicalGroup(2, physical_group):
        _, line = gmsh.model.getAdjacencies(2, tag)
        for i in line:
            _, pts = gmsh.model.getAdjacencies(1, i)
            points.append(pts)
    
    points = list(set([item for sublist in points for item in sublist])) # flatten and filter out any duplicate points in tag list

    # change the formatting to be (dim, tag) list
    entities = []
    for tag in points:
        entities.append((0, tag))
    
    gmsh.model.mesh.setSize(entities, mesh_size)


######## FUSELAGE ########

bottom = gmsh.model.occ.addCone(0, 0, 0, 0, 0, 2 * height / 5, radius - 0.5, radius)
top = gmsh.model.occ.addCylinder(0, 0, 2 * height / 5 + tolerance, 0, 0, 3 * height / 5 + tolerance, radius)


######## TANK ########

tank_hole = gmsh.model.occ.addCylinder(radius - 0.5, 0, 0.4, 0, 0, 2 * height / 5 - 1.3, tank_radius + 0.1)
tank_hole_list = [(3, tank_hole)]

for index in range(0,3):
    tank_hole_list.append(gmsh.model.occ.copy([tank_hole_list[-1]])[0])
    gmsh.model.occ.rotate([tank_hole_list[-1]], 0, 0, 0, 0, 0, 1, math.pi/2)

gmsh.model.occ.cut([(3, bottom)], tank_hole_list)

tank = gmsh.model.occ.addCylinder(radius - 0.5, 0, 0.5, 0, 0, 2 * height / 5 - 1.5, tank_radius)
tank_list = [(3, tank)]
for index in range(0,3):
    tank_list.append(gmsh.model.occ.copy([tank_list[-1]])[0])
    gmsh.model.occ.rotate([tank_list[-1]], 0, 0, 0, 0, 0, 1, math.pi / 2)


######## LANDING LEGS ########

leg = gmsh.model.occ.addCylinder(radius - 0.2, 0, 2 * height / 5 - 0.2, 2.5, 0, - (2 * height / 5 + 1.87), leg_radius)
cone = gmsh.model.occ.addCone(0, 0, 0, 0, 0, 2 * height / 5, radius - 0.5 + tolerance, radius + tolerance)

gmsh.model.occ.cut([(3, leg)], [(3, cone)])

leg_list = [(3, leg)]
gmsh.model.occ.rotate([leg_list[-1]], 0, 0, 0, 0, 0, 1, math.pi / 4)

for index in range(0,3):
    leg_list.append(gmsh.model.occ.copy([leg_list[-1]])[0])
    gmsh.model.occ.rotate([leg_list[-1]], 0, 0, 0, 0, 0, 1, math.pi / 2)


######## BOUNDARY & PHYSICAL GROUPS ########

gmsh.model.occ.synchronize()

# get the surfaces for the two halves of the fuselage and then assign them to physical groups
_ , bottom_surfaces = gmsh.model.getAdjacencies(3, bottom)
_ , top_surfaces = gmsh.model.getAdjacencies(3, top)
ps_top = gmsh.model.addPhysicalGroup(2, top_surfaces, name="Fuselage Top")
ps_bottom = gmsh.model.addPhysicalGroup(2, bottom_surfaces, name="Fuselage Bottom")


# get the surfaces for the tanks and assign them to their own individual physical groups
i = 0
ps_tank_list = []
for volume in tank_list:
    i = i + 1
    up, surface_tags = gmsh.model.getAdjacencies(*volume)
    name = "Tank " + str(i)
    ps_tank_list.append(gmsh.model.addPhysicalGroup(2, surface_tags, name=name))

# get the surfaces for the landing legs and assign them to their own individual physical groups
i = 0
ps_leg_list = []
for volume in leg_list:
    i = i + 1
    up, surface_tags = gmsh.model.getAdjacencies(*volume)
    name = "Leg " + str(i)
    ps_leg_list.append(gmsh.model.addPhysicalGroup(2, surface_tags, name=name))

gmsh.model.occ.synchronize()

# create cylindrical boundary, assign physical group, and then create physical volume
volumes = gmsh.model.occ.getEntities(3)
boundary = gmsh.model.occ.addCylinder(0, 0, -2.2, 0, 0, boundary_radius, boundary_height)
gmsh.model.occ.synchronize()

_ , boundary_surfaces = gmsh.model.getAdjacencies(3, boundary)
gmsh.model.occ.cut([(3, boundary)], volumes)
gmsh.model.occ.synchronize()

ps_ground = gmsh.model.addPhysicalGroup(2, [boundary_surfaces[2]], name="Ground")
ps_space = gmsh.model.addPhysicalGroup(2, [boundary_surfaces[0], boundary_surfaces[1]], name="Space")
pv = gmsh.model.addPhysicalGroup(3, [boundary], name="Volume")

# assign the proper mesh sizes for every physical group.
# NOTE: since cylinders have 3 surfaces but are only defined by 2 points, changing the mesh size for ps_space 
# also changes it for ps_ground. to have a different mesh size for ps_ground it has to be set after the mesh 
# size for ps_space is set.
setMeshSize(ps_bottom, meshsize_lowerfuselage)
setMeshSize(ps_top, meshsize_upperfuselage)
setMeshSize(ps_space, meshsize_space) 
setMeshSize(ps_ground, meshsize_ground)
[setMeshSize(item, meshsize_tanks) for item in ps_tank_list]
[setMeshSize(item, meshsize_landinglegs) for item in ps_leg_list]

gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2) # save msh in ASCII 2 format
gmsh.model.mesh.generate(3)
gmsh.write("blue_moon.msh") # write .msh file 
#gmsh.write("moon.brep") # save .brep file