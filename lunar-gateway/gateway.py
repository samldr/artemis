import gmsh
import math

def setMeshSize(physical_group, mesh_size):
    # function to set the mesh size of surfaces from its physical group tag
    
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


# GLOBAL VARIABLES

tol = 0.01 # the spacing between different physical groups
boundary_radius = 85

docking_radius = 1.3 / 2 # these are used across almost every module so its global
docking_length = 0.17


# MODULE FUNCTIONS

def ppe(a, b, c):
    # References:
    # https://rsdo.gsfc.nasa.gov/images/catalog-rapidIV/Maxar_1300_Data_sheet-Rapid_IV.pdf

    # geometry declarations
    width = 2.4
    depth = 2.2
    height = 3.1

    arm_length = 0.1
    arm_protrusion = 1.5

    panel_width = 2.4
    panel_protrusion = 10


    # mesh declarations
    ms_module = 0.1 * depth
    ms_panel = 0.1 * panel_protrusion


    a = a - width/2
    c = c - depth/2

    # geometry 
    module = gmsh.model.occ.addBox(a, b, c, width, height, depth)

    arm1 = gmsh.model.occ.addBox(a + width/2 - arm_length/2, b + height/2 - arm_length/2, c, arm_length, arm_length, -arm_protrusion)
    arm2 = gmsh.model.occ.copy([(3,arm1)])
    gmsh.model.occ.rotate(arm2, a + width/2, b, c + depth/2, 0, 1, 0, math.pi)

    cyl1 = gmsh.model.occ.addCylinder(a + width/2, b + height, c + depth/2 , 0, docking_length, 0, docking_radius)
    gmsh.model.occ.fuse([(3,module)], [(3, arm1), *arm2, (3, cyl1)])

    panel1 = gmsh.model.occ.addBox(a + width/2 - arm_length/2, b + height/2 - panel_width/2, c - arm_protrusion - tol, arm_length,  panel_width, -panel_protrusion)
    panel2 = gmsh.model.occ.copy([(3,panel1)])
    gmsh.model.occ.rotate(panel2, a + width/2, b, c + depth/2, 0, 1, 0, math.pi)


    # physical grouping
    gmsh.model.occ.synchronize()

    _, module_tags = gmsh.model.getAdjacencies(3, module)
    ps_ppe = gmsh.model.addPhysicalGroup(2, module_tags, name="PPE")
    
    _, panel1_tags = gmsh.model.getAdjacencies(3, panel1)
    ps_ppe_sp1 = gmsh.model.addPhysicalGroup(2, panel1_tags, name="PPE Panel 1")
    
    _, panel2_tags = gmsh.model.getAdjacencies(3, panel2[0][1])
    ps_ppe_sp2 = gmsh.model.addPhysicalGroup(2, panel2_tags, name="PPE Panel 2")


    # mesh sizes
    setMeshSize(ps_ppe, ms_module)
    setMeshSize(ps_ppe_sp1, ms_panel)
    setMeshSize(ps_ppe_sp2, ms_panel)
    
    global dim_ppe  
    dim_ppe = [height + docking_length]

    return [(3, module), (3, panel1), *panel2] # return dimtags

def halo(a, b, c):


    radius = 1.5
    length = 6.1
    
    slope_length = 0.3

    ms_halo = 0.1 * radius

    b = b + docking_length
    module = gmsh.model.occ.addCylinder(a, b + slope_length, c, 0, length - 2 * slope_length, 0, radius)
    cone1 = gmsh.model.occ.addCone(a, b, c, 0, slope_length, 0, docking_radius, radius)
    cone2 = gmsh.model.occ.addCone(a, b + length - slope_length, c, 0, slope_length, 0, radius, docking_radius)
    cyl1 = gmsh.model.occ.addCylinder(a, b, c, 0, -docking_length, 0, docking_radius)
    cyl2 = gmsh.model.occ.addCylinder(a, b + length, c, 0, docking_length, 0, docking_radius)
    cyl3 = gmsh.model.occ.addCylinder(a, b+length/2, c, radius + 2 * docking_length, 0, 0, docking_radius)
    cyl4 = gmsh.model.occ.addCylinder(a, b+length/2, c, -(radius + 2 * docking_length), 0, 0, docking_radius)

    gmsh.model.occ.fuse([(3, module)], [(3, cone1), (3, cone2), (3,cyl1), (3,cyl2), (3,cyl3), (3,cyl4)])

    gmsh.model.occ.synchronize()
    _, module_tags = gmsh.model.getAdjacencies(3, module)
    ps_halo = gmsh.model.addPhysicalGroup(2, module_tags, name="HALO")

    

    setMeshSize(ps_halo, ms_halo)

    global dim_halo
    dim_halo = [length + 2 * docking_length, radius + 2 * docking_length, length / 2 + docking_length]

    return [(3, module)]

def ihab(a, b, c):

    radius = 3.6 / 2 # this is a guess, inner diameter is 3.4m
    length = 6.1

    arm_length = 0.05
    arm_protrusion = 1

    panel_width = 2
    panel_protrusion = 6

    ms_ihab = 0.1 * radius
    ms_panel = 0.1 * panel_protrusion

    b = b + docking_length

    module = gmsh.model.occ.addCylinder(a, b, c, 0, length, 0, radius)
    cyl1 = gmsh.model.occ.addCylinder(a, b, c, 0, -docking_length, 0, docking_radius)
    cyl2 = gmsh.model.occ.addCylinder(a, b + length, c, 0, docking_length, 0, docking_radius)
    cyl3 = gmsh.model.occ.addCylinder(a, b+length/2, c, radius + 2 * docking_length, 0, 0, docking_radius)
    cyl4 = gmsh.model.occ.addCylinder(a, b+length/2, c, -(radius + 2 * docking_length), 0, 0, docking_radius)

    arm1 = gmsh.model.occ.addBox(a - arm_length/2, b + 3/4 * length - arm_length/2, c + radius - 0.2, arm_length, arm_length, arm_protrusion + 0.2)
    arm2 = gmsh.model.occ.copy([(3, arm1)])
    gmsh.model.occ.rotate(arm2, a, b, c, 0, 1, 0, math.pi)

    gmsh.model.occ.fuse([(3, module)], [(3,cyl1), (3,cyl2), (3,cyl3), (3,cyl4), (3, arm1), *arm2])

    panel1 = gmsh.model.occ.addBox(a - arm_length/2, b + 3/4*length - panel_width/2, c + tol + radius + arm_protrusion, arm_length,  panel_width, panel_protrusion)
    panel2 = gmsh.model.occ.copy([(3, panel1)])
    gmsh.model.occ.rotate(panel2, a, b, c, 0, 1, 0, math.pi)


    gmsh.model.occ.synchronize()
    _, module_tags = gmsh.model.getAdjacencies(3, module)
    ps_ihab = gmsh.model.addPhysicalGroup(2, module_tags, name="I-HAB")

    _, panel1_tags = gmsh.model.getAdjacencies(3, panel1)
    ps_ihab_sp1 = gmsh.model.addPhysicalGroup(2, panel1_tags, name="I-HAB Panel 1")
    
    _, panel2_tags = gmsh.model.getAdjacencies(3, panel2[0][1])
    ps_ihab_sp2 = gmsh.model.addPhysicalGroup(2, panel2_tags, name="I-HAB Panel 2")

    setMeshSize(ps_ihab, ms_ihab)
    setMeshSize(ps_ihab_sp1, ms_panel)
    setMeshSize(ps_ihab_sp2, ms_panel)

    global dim_ihab
    dim_ihab = [length + 2 * docking_length, radius + 2 * docking_length, length / 2 + docking_length]

    return [(3, module), (3, panel1), *panel2] # return dimtags

def orion(a, b, c):
    # find a better source for these measurements

    crew_length = 3.3528
    crew_radius = 5.0292 / 2
    
    service_length = crew_length #4.78536 
    service_radius =  4 / 2 # this is a guess based on the dimensions in https://www.esa.int/Science_Exploration/Human_and_Robotic_Exploration/Orion/Artemis_IV

    heatshield_thickness = 1

    arm_length = 0.05
    arm_protrusion = 1

    panel_protrusion = 7  - arm_protrusion
    panel_width = 2

    ms_orion = 0.1 * crew_radius
    ms_panel = 0.1 * panel_protrusion

    module = gmsh.model.occ.addCylinder(a, b, c, 0, 2 * docking_length, 0, docking_radius)
    cyl1 = gmsh.model.occ.addCone(a, b + 2 * docking_length, c, 0, crew_length, 0, docking_radius, crew_radius)
    cyl2 = gmsh.model.occ.addCylinder(a, b + 2 * docking_length + crew_length, c, 0, service_length, 0, service_radius)
    cyl3 = gmsh.model.occ.addCylinder(a, b + 2 * docking_length + crew_length - heatshield_thickness/2, c, 0, heatshield_thickness, 0, crew_radius)
    cyl4 = gmsh.model.occ.addCone(a, b + 2 * docking_length + crew_length + service_length, c, 0, 1, 0, 0.5, 1)

    arm1 = gmsh.model.occ.addBox(a - arm_length/2, b + crew_length + 7/8 * service_length - arm_length/2, c + service_radius - 0.2, arm_length, arm_length, arm_protrusion + 0.2)
    gmsh.model.occ.rotate([(3, arm1)], a, b, c, 0, 1, 0, math.pi/3)

    arm2 = gmsh.model.occ.copy([(3, arm1)])
    gmsh.model.occ.rotate(arm2, a, b, c, 0, 1, 0, math.pi)

    arm3 = gmsh.model.occ.copy([(3, arm1)])
    gmsh.model.occ.rotate(arm3, a, b, c, 0, 1, 0, 4 * math.pi/3)
    
    arm4 = gmsh.model.occ.copy(arm3)
    gmsh.model.occ.rotate(arm4, a, b, c, 0, 1, 0, math.pi)

    gmsh.model.occ.fuse([(3, module)], [(3,cyl1), (3,cyl2), (3,cyl3), (3,cyl4), (3, arm1), *arm2, *arm3, *arm4])

    panel1 = gmsh.model.occ.addBox(a - arm_length/2, b + crew_length + 7/8 * service_length - panel_width/2, c + service_radius + tol + arm_protrusion, arm_length,  panel_width, panel_protrusion)
    gmsh.model.occ.rotate([(3, panel1)], a, b + crew_length + 7/8 * service_length, c + service_radius + tol + arm_protrusion, 0, 0, 1, math.pi/2)
    gmsh.model.occ.rotate([(3, panel1)], a, b, c, 0, 1, 0, math.pi/3)

    panel2 = gmsh.model.occ.copy([(3, panel1)])
    gmsh.model.occ.rotate(panel2, a, b, c, 0, 1, 0, math.pi)

    panel3 = gmsh.model.occ.copy([(3, panel1)])
    gmsh.model.occ.rotate(panel3, a, b, c, 0, 1, 0, 4 * math.pi/3)
    
    panel4 = gmsh.model.occ.copy(panel3)
    gmsh.model.occ.rotate(panel4, a, b, c, 0, 1, 0, math.pi)

    gmsh.model.occ.synchronize()
    _, module_tags = gmsh.model.getAdjacencies(3, module)
    ps_orion = gmsh.model.addPhysicalGroup(2, module_tags, name="Orion")

    _, panel1_tags = gmsh.model.getAdjacencies(3, panel1)
    ps_orion_sp1 = gmsh.model.addPhysicalGroup(2, panel1_tags, name="Orion Panel 1")
    
    _, panel2_tags = gmsh.model.getAdjacencies(3, panel2[0][1])
    ps_orion_sp2 = gmsh.model.addPhysicalGroup(2, panel2_tags, name="Orion Panel 2")
    
    _, panel3_tags = gmsh.model.getAdjacencies(3, panel3[0][1])
    ps_orion_sp3 = gmsh.model.addPhysicalGroup(2, panel3_tags, name="Orion Panel 3")

    _, panel4_tags = gmsh.model.getAdjacencies(3, panel4[0][1])
    ps_orion_sp4 = gmsh.model.addPhysicalGroup(2, panel4_tags, name="Orion Panel 4")

    setMeshSize(ps_orion, ms_orion)
    setMeshSize(ps_orion_sp1, ms_panel)
    setMeshSize(ps_orion_sp2, ms_panel)
    setMeshSize(ps_orion_sp3, ms_panel)
    setMeshSize(ps_orion_sp4, ms_panel)

    global dim_orion
    dim_orion = [2 * docking_length + crew_length + service_length]

    return [(3, module), (3, panel1), *panel2, *panel3, *panel4] # return dimtags

def esprit(a, b, c):
    length = 6.4 
    radius = 4.6 / 2
    hex_length = 2.5
    smaller_radius = math.sqrt(3)/2 * radius

    ms_esprit = 0.1 * smaller_radius

    module = gmsh.model.occ.addCylinder(a, b, c, 2* docking_length, 0, 0, docking_radius)

    points = [(0, gmsh.model.occ.addPoint(a + 2*docking_length, b, c + radius))]
    
    for i in range(1, 6):
        points.append(gmsh.model.occ.copy([points[-1]])[0])
        gmsh.model.occ.rotate([points[-1]], a, b, c, 1, 0, 0, math.pi/3)

    lines = []
    for i in range(1, 6):
        lines.append(gmsh.model.occ.addLine(points[i - 1][1], points[i][1]))
    lines.append(gmsh.model.occ.addLine(points[-1][1], points[0][1]))
    cl = gmsh.model.occ.addCurveLoop(lines)

    s = gmsh.model.occ.addPlaneSurface([cl])
    hex = gmsh.model.occ.extrude([(2, s)], hex_length, 0, 0)[1]

    cyl1 = gmsh.model.occ.addCylinder(a + 2 * docking_length + hex_length, b, c, length - hex_length - 3 * docking_length, 0, 0, smaller_radius)
    cyl2 = gmsh.model.occ.addCylinder(a + length - docking_length, b, c, docking_length, 0, 0, docking_radius)
    
    gmsh.model.occ.fuse([(3, module)], [hex, (3, cyl1), (3, cyl2)])

    gmsh.model.occ.rotate([(3, module)], a, b, c, 0, 0, 1, math.pi)

    gmsh.model.occ.synchronize()
    _, module_tags = gmsh.model.getAdjacencies(3, module)
    ps_esprit = gmsh.model.addPhysicalGroup(2, module_tags, name="ESPRIT")

    setMeshSize(ps_esprit, ms_esprit)

    global dim_esprit
    dim_esprit = [length + docking_length]

    return [(3, module)]

def bluemoon(a, b, c):
    # this part of the geometry is taken from the blue_moon.py script

    ######## MODEL PARAMETERS ########

    height = 16 # overall lander height
    radius = 3 # fuselage radius
    tank_radius = 1.25 # radial tank radius

    c = c - height - docking_length - 2 * tol

    ######## MESH PARAMETERS ########

    meshsize_lowerfuselage = 0.1 * radius # upper fuselage
    meshsize_upperfuselage = 0.1 * (radius - 0.5) # lower fuselage
    meshsize_tanks = 0.1 * tank_radius # radially mounted tanks

    ######## FUSELAGE ########

    bottom = gmsh.model.occ.addCone(a, b, c, 0, 0, 2 * height / 5, radius - 0.5, radius)
    top = gmsh.model.occ.addCylinder(a, b, c+ 2 * height / 5 + tol, 0, 0, 3 * height / 5 + tol, radius)

    cyl1 = gmsh.model.occ.addCylinder(a, b, c + height + 2 * tol, 0, 0, docking_length, docking_radius)

    gmsh.model.occ.fuse([(3, top)], [(3, cyl1)])

    ######## TANK ########

    tank_hole = gmsh.model.occ.addCylinder(a + radius - 0.5, b, c + 0.4, 0, 0, 2 * height / 5 - 1.3, tank_radius + 0.1)
    tank_hole_list = [(3, tank_hole)]

    for index in range(0,3):
        tank_hole_list.append(gmsh.model.occ.copy([tank_hole_list[-1]])[0])
        gmsh.model.occ.rotate([tank_hole_list[-1]], a, b, c, 0, 0, 1, math.pi/2)

    gmsh.model.occ.cut([(3, bottom)], tank_hole_list)

    tank = gmsh.model.occ.addCylinder(a + radius - 0.5, b , c + 0.5, 0, 0, 2 * height / 5 - 1.5, tank_radius)
    tank_list = [(3, tank)]
    for index in range(0,3):
        tank_list.append(gmsh.model.occ.copy([tank_list[-1]])[0])
        gmsh.model.occ.rotate([tank_list[-1]], a, b, c, 0, 0, 1, math.pi / 2)



    gmsh.model.occ.rotate([*tank_list, (3, bottom), (3, top)], a, b, c + height + docking_length + 2 * tol, 0, 1, 0, -math.pi/2)

    ######## BOUNDARY & PHYSICAL GROUPS ########

    gmsh.model.occ.synchronize()

    # get the surfaces for the two halves of the fuselage and then assign them to physical groups
    _ , bottom_surfaces = gmsh.model.getAdjacencies(3, bottom)
    _ , top_surfaces = gmsh.model.getAdjacencies(3, top)
    ps_top = gmsh.model.addPhysicalGroup(2, top_surfaces, name="Blue Moon Fuselage Top")
    ps_bottom = gmsh.model.addPhysicalGroup(2, bottom_surfaces, name="Blue Moon Fuselage Bottom")

    # get the surfaces for the tanks and assign them to their own individual physical groups
    i = 0
    ps_tank_list = []
    for volume in tank_list:
        i = i + 1
        _, surface_tags = gmsh.model.getAdjacencies(*volume)
        name = "Blue Moon Tank " + str(i)
        ps_tank_list.append(gmsh.model.addPhysicalGroup(2, surface_tags, name=name))

    setMeshSize(ps_bottom, meshsize_lowerfuselage)
    setMeshSize(ps_top, meshsize_upperfuselage)
    [setMeshSize(item, meshsize_tanks) for item in ps_tank_list]
    
    return [*tank_list, (3, bottom), (3, top)]

def dragonxl(a, b, c):
    radius = 1.5
    length = 6.1
    
    slope_length = 0.3

    back_length = 1.5
    back_radius = 1.2

    arm_length = 0.1
    arm_protrusion = 0.75

    panel_width = 2
    panel_protrusion = 7.5

    ms_dragonxl = 0.1 * radius
    ms_panel = 0.1 * panel_protrusion

    module = gmsh.model.occ.addCylinder(a, b + slope_length, c, 0, length - 2 * slope_length - back_length, 0, radius)
    cyl1 = gmsh.model.occ.addCylinder(a, b, c, 0, -docking_length, 0, docking_radius)
    cyl2 = gmsh.model.occ.addCylinder(a, b+ length - 2 * slope_length - back_length, c, 0, back_length, 0, back_radius)
    cone1 = gmsh.model.occ.addCone(a, b, c, 0, slope_length, 0, docking_radius, radius)


    arm1 = gmsh.model.occ.addBox(a - arm_length/2, b + length - back_length + 1/2 * back_length - arm_length/2, c + back_radius - 0.2, arm_length, arm_length, arm_protrusion + 0.2)
    arm2 = gmsh.model.occ.copy([(3, arm1)])
    gmsh.model.occ.rotate(arm2, a, b, c, 0, 1, 0, math.pi)

    gmsh.model.occ.fuse([(3, module)], [(3, cone1), (3, cyl1), (3, cyl2), (3, arm1), *arm2])

    panel1 = gmsh.model.occ.addBox(a - panel_width/2, b + length - back_length + 1/2 * back_length - arm_length/2, c + tol + back_radius + arm_protrusion, panel_width, arm_length, panel_protrusion)
    panel2 = gmsh.model.occ.copy([(3, panel1)])
    gmsh.model.occ.rotate(panel2, a, b, c, 0, 1, 0, math.pi)
    

    gmsh.model.occ.rotate([(3, module), (3, panel1), *panel2], a, b, c, 0, 0, 1, math.pi/2)

    gmsh.model.occ.synchronize()
    _, module_tags = gmsh.model.getAdjacencies(3, module)
    ps_dragonxl = gmsh.model.addPhysicalGroup(2, module_tags, name="Dragon XL")

    _, panel1_tags = gmsh.model.getAdjacencies(3, panel1)
    ps_dragonxl_sp1 = gmsh.model.addPhysicalGroup(2, panel1_tags, name="Dragon XL Panel 1")
    
    _, panel2_tags = gmsh.model.getAdjacencies(3, panel2[0][1])
    ps_dragonxl_sp2 = gmsh.model.addPhysicalGroup(2, panel2_tags, name="Dragon XL Panel 2")

    

    setMeshSize(ps_dragonxl, ms_dragonxl)
    setMeshSize(ps_dragonxl_sp1, ms_panel)
    setMeshSize(ps_dragonxl_sp2, ms_panel)


    return [(3, module), (3, panel1), *panel2]

def airlock(a, b, c):
    radius  = 2.5/2
    length = 3.5 

    ms_airlock = 0.1 * radius
    
    module = gmsh.model.occ.addCylinder(a, b, c, 0, docking_length, 0, docking_radius)
    cyl1 = gmsh.model.occ.addCylinder(a, b + docking_length, c, 0, length, 0, radius)

    gmsh.model.occ.fuse([(3, module)], [(3, cyl1)])

    gmsh.model.occ.rotate([(3, module)], a, b, c, 0, 0, 1, math.pi/2)

    gmsh.model.occ.synchronize()
    _, module_tags = gmsh.model.getAdjacencies(3, module)
    ps_airlock = gmsh.model.addPhysicalGroup(2, module_tags, name="Airlock")

    setMeshSize(ps_airlock, ms_airlock)

    return [(3, module)]


# CREATE GEOMETRY

gmsh.initialize()

offset = -11.6128 

ppe_volumes = ppe(0, offset, 0)
halo_volumes = halo(0, offset + dim_ppe[0] + tol, 0)
ihab_volumes = ihab(0, offset + dim_ppe[0] + dim_halo[0] + 2 * tol, 0)
orion_volumes = orion(0, offset + dim_ppe[0] + dim_halo[0] + dim_ihab[0] + 3 * tol, 0)
bluemoon_volumes = bluemoon(dim_halo[1] + tol, offset + dim_halo[2] + dim_ppe[0] + tol, 0)
esprit_volumes = esprit(-(dim_halo[1] + tol), offset + dim_halo[2] + dim_ppe[0] + tol, 0)
dragonxl_volumes = dragonxl(-(dim_halo[1] + dim_esprit[0] + 2 * tol ), offset + dim_halo[2] + dim_ppe[0] + tol, 0)
airlock_volumes = airlock(-(dim_ihab[1] + tol), offset + dim_ppe[0] + dim_halo[0] + dim_ihab[2] + 2 * tol, 0)

station_volumes = [*ppe_volumes, *halo_volumes, *ihab_volumes, *orion_volumes, *bluemoon_volumes, *esprit_volumes, *dragonxl_volumes, *airlock_volumes]

# the offset to center the station. all the modules need to be created to find the length so you must run the script with the 2 lines below uncommented
# to find the length, and then take the printed value and replace the offset declaration at the top with it.

# offset = -(dim_ppe[0] + dim_halo[0] + dim_ihab[0] + dim_orion[0] + 3 * tol)/2
# print(offset)

boundary = gmsh.model.occ.addSphere(0, 0, 0, boundary_radius)

gmsh.model.occ.synchronize()
_ , boundary_surfaces = gmsh.model.getAdjacencies(3, boundary)
gmsh.model.occ.cut([(3, boundary)], station_volumes)
gmsh.model.occ.synchronize()

ps_space = gmsh.model.addPhysicalGroup(2, [*boundary_surfaces], name="Space")
pv = gmsh.model.addPhysicalGroup(3, [boundary], name="Volume")


# MESHING

setMeshSize(ps_space, 0.1 * boundary_radius) 

gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2) # save msh in ASCII 2 format
gmsh.write("gateway.brep")
gmsh.model.mesh.generate(3)
gmsh.write("gateway.msh")