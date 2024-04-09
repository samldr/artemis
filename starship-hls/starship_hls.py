import gmsh
from math import pi

gmsh.initialize()


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


######## MODEL PARAMETERS ########

# fuselage
fuselage_radius = 4.5
fuselage_height = 36.28

# nosecone
nosecone_height = 13.72 
nc1 = 0.75 # parameter for nosecone bezier curve
nc2 = 2 # parameter for nosecone bezier curve

# engines
engine_bay_height = 7.32
engine_bay_thickness = 0.5
engine_radius = 0.65

# landing legs
landing_leg_housing_width = 2
landing_leg_housing_height = 10
landing_leg_length = 8

# bounding box
boundary_height = 80
boundary_radius = 30

# spacing for different physical groups
spacing = 0.1

######## MESH PARAMETERS ########


meshsize_fuselage = 0.1 * fuselage_radius
meshsize_solarpanels = 0.2
meshsize_landinglegs = 0.2
meshsize_lunarsurface = 0.1 * boundary_radius
meshsize_space = 0.1 * boundary_radius



# meshsize_fuselage = 0.2 * fuselage_radius
# meshsize_solarpanels = 1
# meshsize_landinglegs = 0.025
# meshsize_lunarsurface = 3
# meshsize_space = 3



######## FUSELAGE ########

lander = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, fuselage_height, fuselage_radius)

# nose cone
p1 = gmsh.model.occ.addPoint(fuselage_radius, 0, fuselage_height)
p2 = gmsh.model.occ.addPoint(fuselage_radius * nc1, 0, fuselage_height + nosecone_height)
p3 = gmsh.model.occ.addPoint(0, 0, fuselage_height + nosecone_height)
c1 = gmsh.model.occ.addBezier([p1, p2, p3])

s1 = gmsh.model.occ.revolve([(1, c1)], 0, 0, 0, 0, 0, 1, 2*pi)
v1 = gmsh.model.occ.extrude(s1, 0, 0, -nosecone_height)


nosecone = None
for dimtag in v1:
    if dimtag[0] == 3:
        nosecone = dimtag
        break
gmsh.model.occ.fuse([(3, lander)], [nosecone])
gmsh.model.occ.remove(v1, recursive=True)

gmsh.model.occ.synchronize()
gmsh.write("fuselage.brep")

# engine bay
v1 = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, engine_bay_height, fuselage_radius - engine_bay_thickness) # empty space for engine bay
v2 = gmsh.model.occ.addCone(0, 0, engine_bay_height, 0, 0, - engine_bay_height / 2, fuselage_radius - engine_bay_thickness, 0) # add small inverted cone at the top of the engine bay
gmsh.model.occ.cut([(3,lander)], [(3, v1)]) # remove the space from the fuselage
gmsh.model.occ.fuse([(3,lander)], [(3, v2)]) # fuse the cone with the fuselage


######## ENGINES ########

large_engine = gmsh.model.occ.addCone(2/3 * fuselage_radius - engine_bay_thickness, 0, 0, 0, 0, engine_bay_height, 2 * engine_radius, engine_radius)
small_engine = gmsh.model.occ.addCone(1/3 * fuselage_radius - engine_bay_thickness, 0, 0, 0, 0, engine_bay_height, engine_radius, 1/2 * engine_radius)

large_engine_list = [(3, large_engine)]
small_engine_list = [(3, small_engine)]

gmsh.model.occ.rotate([small_engine_list[-1]], 0, 0, 0, 0, 0, 1, pi/3)

for index in range(0, 2):
    small_engine_list.append(gmsh.model.occ.copy([small_engine_list[-1]])[0])
    gmsh.model.occ.rotate([small_engine_list[-1]], 0, 0, 0, 0, 0, 1, 2*pi/3)
    large_engine_list.append(gmsh.model.occ.copy([large_engine_list[-1]])[0])
    gmsh.model.occ.rotate([large_engine_list[-1]], 0, 0, 0, 0, 0, 1, 2*pi/3)


temp =  gmsh.model.occ.fuse([(3,lander)], small_engine_list)
lander = temp[0][0][1]

gmsh.model.occ.fuse([(3,lander)], large_engine_list)

######## SOLAR PANELS ########

v1 = gmsh.model.occ.addCylinder(0, 0, fuselage_height*13/16 , 0, 0, 6, fuselage_radius + 0.2, angle = 0.8*pi/4)
v2 = gmsh.model.occ.addCylinder(0, 0, fuselage_height*13/16 , 0, 0, 6, fuselage_radius + 0.1, angle = 0.8*pi/4)
gmsh.model.occ.cut([(3, v1)], [(3, v2)])

solar_panel_list = [(3,v1)]

for index in range(0, 7):
    solar_panel_list.append(gmsh.model.occ.copy([solar_panel_list[-1]])[0])
    gmsh.model.occ.rotate([solar_panel_list[-1]], 0, 0, 0, 0, 0, 1, 2*pi/8)


######## LANDING LEG HOUSING ########

# points
p1 = gmsh.model.occ.addPoint(landing_leg_housing_width, fuselage_radius - engine_bay_thickness, 0)
p2 = gmsh.model.occ.addPoint(-landing_leg_housing_width, fuselage_radius - engine_bay_thickness, 0)
p3 = gmsh.model.occ.addPoint(-landing_leg_housing_width, fuselage_radius - engine_bay_thickness, landing_leg_housing_height)
p4 = gmsh.model.occ.addPoint(landing_leg_housing_width, fuselage_radius - engine_bay_thickness, landing_leg_housing_height)

p5 = gmsh.model.occ.addPoint(landing_leg_housing_width - 0.2, fuselage_radius + engine_bay_thickness, 0)
p6 = gmsh.model.occ.addPoint(-(landing_leg_housing_width - 0.2), fuselage_radius + engine_bay_thickness, 0)
p7 = gmsh.model.occ.addPoint(-(landing_leg_housing_width - 0.9), fuselage_radius + engine_bay_thickness, 3/5 * landing_leg_housing_height)
p8 = gmsh.model.occ.addPoint(landing_leg_housing_width - 0.9, fuselage_radius + engine_bay_thickness, 3/5 * landing_leg_housing_height)

# lines
l1 = gmsh.model.occ.addLine(p1, p2)
l2 = gmsh.model.occ.addLine(p2, p3)
l3 = gmsh.model.occ.addLine(p3, p4)
l4 = gmsh.model.occ.addLine(p4, p1)

l5 = gmsh.model.occ.addLine(p5, p6)
l6 = gmsh.model.occ.addLine(p6, p7)
l7 = gmsh.model.occ.addLine(p7, p8)
l8 = gmsh.model.occ.addLine(p8, p5)

l9 = gmsh.model.occ.addLine(p1, p5)
l10 = gmsh.model.occ.addLine(p2, p6)
l11 = gmsh.model.occ.addLine(p3, p7)
l12 = gmsh.model.occ.addLine(p4, p8)

# curve loops
c1 = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
c2 = gmsh.model.occ.addCurveLoop([l5, l6, l7, l8])
c3 = gmsh.model.occ.addCurveLoop([l9, l5, -l10, -l1])
c4 = gmsh.model.occ.addCurveLoop([l11, l7, -l12, -l3])
c5 = gmsh.model.occ.addCurveLoop([l4, l9, -l8, -l12])
c6 = gmsh.model.occ.addCurveLoop([l2, l11, -l6, -l10])

# surfaces
s1 = gmsh.model.occ.addPlaneSurface([c1])
s2 = gmsh.model.occ.addPlaneSurface([c2])
s3 = gmsh.model.occ.addPlaneSurface([c3])
s4 = gmsh.model.occ.addPlaneSurface([c4])
s5 = gmsh.model.occ.addSurfaceFilling(c5)
s6 = gmsh.model.occ.addSurfaceFilling(c6)

# volume
sl1 = gmsh.model.occ.addSurfaceLoop([s1, s3, s2, s4, s5, s6])

# add the indent using extrude and cut
v1 = gmsh.model.occ.addVolume([sl1])
v2 = gmsh.model.occ.extrude([(2, s2)], 0, -0.5, 0)

indent = None
for dimtag in v2:
    if dimtag[0] == 3:
        indent = dimtag
        break

gmsh.model.occ.cut([(3, v1)],[indent])[0]
gmsh.model.occ.remove(v2, recursive=True)

housing_list = [(3,v1)]

for index in range(1,4):
    housing_list.append(gmsh.model.occ.copy([housing_list[-1]])[0])
    gmsh.model.occ.rotate([housing_list[-1]], 0, 0, 0, 0, 0, 1, pi/2)

gmsh.model.occ.fuse([(3,lander)], housing_list)


######## LANDING LEGS ########

landing_gear = gmsh.model.occ.addBox(0 - 0.2, fuselage_radius - 0.5, 2/5 * landing_leg_housing_height + 0.2, 0 + 0.4, landing_leg_length, 0.4)
gmsh.model.occ.rotate([(3, landing_gear)], 0, fuselage_radius, 2/5 * landing_leg_housing_height, 1, 0, 0, -pi/3)
gmsh.model.occ.translate([(3, landing_gear)],0, -0.25, -1/4 * landing_leg_housing_height)


# landing_gear = gmsh.model.occ.addCylinder(0, fuselage_radius - 0.25 + 4, 2/5 * landing_leg_housing_height, 0, fuselage_radius,  -landing_leg_length, 1)

box = gmsh.model.occ.addBox(0 - fuselage_radius - spacing, 0 - fuselage_radius - spacing, -10, 2* fuselage_radius + 2 * spacing, 2 * fuselage_radius + 2 * spacing, 50)

gmsh.model.occ.cut([(3, landing_gear)], [(3, box)])
p1 = gmsh.model.occ.addPoint(landing_leg_housing_width, fuselage_radius - engine_bay_thickness, 0)
p2 = gmsh.model.occ.addPoint(-landing_leg_housing_width, fuselage_radius - engine_bay_thickness, 0)

gmsh.model.occ.synchronize()
pts = gmsh.model.getBoundary([(3,landing_gear)], recursive=True)

pts_coords = []
#pts_coords.append(gmsh.model.getValue(0, pts[0][1], [])) # had to look in gui to find these points
pts_coords.append(gmsh.model.getValue(0, pts[1][1], []))
pts_coords.append(gmsh.model.getValue(0, p1, []))
pts_coords.append(gmsh.model.getValue(0, p2, []))


gmsh.model.occ.synchronize()

# v2 = gmsh.model.occ.addCylinder(pts_coords[0][0], pts_coords[0][1], pts_coords[0][2], -pts_coords[0][0] + landing_leg_housing_width -0.2, -pts_coords[0][1] + (fuselage_radius - engine_bay_thickness) - 0.2, -pts_coords[0][2] + 0.1, 0.1)
# v3 = gmsh.model.occ.addCylinder(pts_coords[0][0], pts_coords[0][1], pts_coords[0][2], -pts_coords[0][0] - landing_leg_housing_width +0.2, -pts_coords[0][1] + (fuselage_radius - engine_bay_thickness) - 0.2, -pts_coords[0][2] + 0.1, 0.1)

p1 = gmsh.model.occ.addPoint(pts_coords[0][0] + 1.25, pts_coords[0][1] + 2, pts_coords[0][2] - 0.25)
p2 = gmsh.model.occ.addPoint(pts_coords[0][0] + 1.5, pts_coords[0][1] - 1, pts_coords[0][2] - 0.25)
p3 = gmsh.model.occ.addPoint(pts_coords[0][0] - 1.25, pts_coords[0][1] + 2, pts_coords[0][2] - 0.25)
p4 = gmsh.model.occ.addPoint(pts_coords[0][0] - 1.5, pts_coords[0][1] - 1, pts_coords[0][2] - 0.25)

l1 = gmsh.model.occ.addLine(p1, p3)
l2 = gmsh.model.occ.addLine(p3, p4)
l3 = gmsh.model.occ.addLine(p4, p2)
l4 = gmsh.model.occ.addLine(p2, p1)

c1 = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
s1 = gmsh.model.occ.addPlaneSurface([c1])

v4 = gmsh.model.occ.extrude([(2, s1)], 0, 0, 0.25)

foot = None
for dimtag in v4:
    if dimtag[0] == 3:
        foot = dimtag
        break

gmsh.model.occ.fuse([(3, landing_gear)], [foot])
gmsh.model.occ.remove(v4, recursive=True)

landing_gear_list = [(3,landing_gear)]

for index in range(1,4):
    landing_gear_list.append(gmsh.model.occ.copy([landing_gear_list[-1]])[0])
    gmsh.model.occ.rotate([landing_gear_list[-1]], 0, 0, 0, 0, 0, 1, pi/2)

######## BOUNDARY & PHYSICAL GROUPS ########

gmsh.model.occ.synchronize()


# lander = 1

# _, lander_tags = gmsh.model.getAdjacencies(3, lander)
# lander_tags = lander_tags + [18, 19, 20, 21, 22]
# ps_lander = gmsh.model.addPhysicalGroup(2, lander_tags, name="Lander")

sp_surfaces = []
ps_solar_panel_list = []
i = 0
for volume in solar_panel_list:
    i = i + 1
    up, surface_tags = gmsh.model.getAdjacencies(*volume)
    sp_surfaces.append(surface_tags)
    name = "Solar Panel " + str(i)
    ps_solar_panel_list.append(gmsh.model.addPhysicalGroup(2, surface_tags, name=name))
    
lg_surfaces = []
ps_landing_gear_list = []   
i = 0    
for volume in landing_gear_list:
    i = i + 1 
    up, surface_tags = gmsh.model.getAdjacencies(*volume)
    lg_surfaces.append(surface_tags)
    name = "Landing Leg " + str(i)
    ps_landing_gear_list.append(gmsh.model.addPhysicalGroup(2, surface_tags, name=name))



volumes = gmsh.model.occ.getEntities(3)
boundary = gmsh.model.occ.addCylinder(0, 0, -5, 0, 0, boundary_height, boundary_radius)
gmsh.model.occ.synchronize()
gmsh.model.occ.synchronize()


_ , boundary_surfaces = gmsh.model.getAdjacencies(3, boundary)
gmsh.model.occ.cut([(3, boundary)], volumes)
gmsh.model.occ.synchronize()


ps_space = gmsh.model.addPhysicalGroup(2, [boundary_surfaces[0],boundary_surfaces[1]], name="Space")
ps_lunar_surface = gmsh.model.addPhysicalGroup(2, [boundary_surfaces[2]], name="Lunar Surface")

all_surfaces_dimtag = gmsh.model.occ.getEntities(2)
all_surfaces = [item[1] for item in all_surfaces_dimtag]


not_lander_surfaces = [item for sublist in [boundary_surfaces, *sp_surfaces, *lg_surfaces] for item in sublist]
lander_surfaces = list(filter(lambda x: x not in not_lander_surfaces, all_surfaces))


ps_lander = gmsh.model.addPhysicalGroup(2, lander_surfaces, name="Lander")

pv = gmsh.model.addPhysicalGroup(3, [boundary], name="Volume")


######## MESHING ########

gmsh.model.occ.synchronize()

setMeshSize(ps_lander, meshsize_fuselage)
[setMeshSize(item, meshsize_solarpanels) for item in ps_solar_panel_list]
[setMeshSize(item, meshsize_landinglegs) for item in ps_landing_gear_list]
setMeshSize(ps_lunar_surface, meshsize_lunarsurface)
setMeshSize(ps_space, meshsize_space)

gmsh.write("starship_hls.brep")

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2) # save msh in ASCII 2 format
gmsh.model.mesh.generate(3)
gmsh.write("starship_hls.msh")
gmsh.finalize()