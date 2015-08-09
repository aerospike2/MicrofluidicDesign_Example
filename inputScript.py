# ------------------------------------------------------------------------- #
# ---------------------- MICROFLUIDIC CHANNEL DESIGN ---------------------- # 
# ------------------------------------------------------------------------- #

# Abhishek Sharma, BioMIP, Ruhr University, Bochum, Germany
# Last Updated : 13th May, 2013


""" 
This script creates microfluidc geometry using SALOME Platform for CFD simulations using OpenFOAM.


                         Description of Geometry 
                        -------------------------
                               0__1
                               |  |
                               |  |
                               |  |
                               |  | 
          12______10           |  |
           |_____   \          |  |                               
          18   16\__ \__       |  |        
                    \__ \__    |  |     
                      14   8___3  2
                        |   ___   |
                     _15 __9   4  5
                  __/ __/      |  |
          19___17/   /         |  |
          |________/ 	       |  |
    	  13      11           |  |
                               |  |
                               |  |
                               |  |
                               |  |
                               7__6
"""

import numpy as np
import math as m
import salome

# --------------------------------------------------------------------------- # 
# Input Parameters

fac = 0.001							# multiplication factor 
widM = 0.1							# width of main channel
lenMT = 1.0							# length of main channel top
widC = 0.04							# width of thin connecting channel
lenMB = 1.0							# length of the main channel bottom
lenC = 0.1							# length of thin connecting channel
lInc1 = 0.5							# length of inclined channel (8-10, 9-11)
lenH = 0.5							# length of horizontal channel (10-12, 13-11)
lenX = 0.1							# length between points (8-14, 9-15) 
delX = 0.1							# change in the length at inclined points
ang = 45.*m.pi/180.						# angle of inclination
hgh = 0.05							# height of the channels
hghS = 0.005							# height of small channel

# --------------------------------------------------------------------------- #
# Creating Points (Two dimensional Projection)

p = []

p0 = [0.0, 0.0]
p.append(p0)

p1 = [p0[0] + widM, p0[1]]
p.append(p1)

p2 = [p1[0], p1[1] - lenMT]  
p.append(p2) 

p3 = [p0[0], p0[1] - lenMT]  
p.append(p3)

p4 = [p3[0], p3[1] - widC]  
p.append(p4)

p5 = [p2[0], p2[1] - widC]  
p.append(p5)

p6 = [p5[0], p5[1] - lenMB]
p.append(p6)

p7 = [p4[0], p4[1] - lenMB]
p.append(p7)

p8 = [p3[0] - lenC, p3[1]]  
p.append(p8) 

p9 = [p4[0] - lenC, p4[1]]  
p.append(p9)

p10 = [p8[0] - lInc1*m.cos(ang), p8[1] + lInc1*m.sin(ang)]  
p.append(p10)

p11 = [p9[0] - lInc1*m.cos(ang), p9[1] - lInc1*m.sin(ang)]
p.append(p11)

p12 = [p10[0] - lenH, p10[1]]  
p.append(p12)

p13 = [p11[0] - lenH, p11[1]]
p.append(p13)

p14 = [p8[0] - lenX, p8[1]]  
p.append(p14)

p15 = [p9[0] - lenX, p9[1]]
p.append(p15)

p16 = [p14[0] - (lInc1 - delX)*m.cos(ang), p14[1] + (lInc1 - delX)*m.sin(ang)]  
p.append(p16)

p17 = [p15[0] - (lInc1 - delX)*m.cos(ang), p15[1] - (lInc1 - delX)*m.sin(ang)]
p.append(p17)

p18 = [p12[0], p16[1]]  
p.append(p18)

p19 = [p13[0], p17[1]]
p.append(p19)

# Printing out the created points
for i in range(len(p)):
   print '('+ repr(p[i][0]) + '  ' + repr(p[i][1]) + '  ' + repr(0.0) + ')      //' + repr(i)
      
for i in range(len(p)):
   print '('+repr(p[i][0]) + '  ' + repr(p[i][1]) + '  ' + repr(hgh) + ')       //' + repr(len(p) + i)
                    
nop = len(p)
arr = np.array(p)
         
# -------------------------------------------------------------------------- #
# Defining regions for setFieldsDict (OpenFOAM)

print repr(fac*arr[13]) + ' '+ repr(fac*arr[10])

# -------------------------------------------------------------------------- #
# Creating geometry in SALOME

import salome
import geompy
gg = salome.ImportComponentGUI("GEOM")

pX = []

for i in range(len(p)):
   pX.append(geompy.MakeVertex(fac*arr[i,0], fac*arr[i,1], fac*0.0))
         
for i in range(len(p)):
   pX.append(geompy.MakeVertex(fac*arr[i,0], fac*arr[i,1], fac*hgh))

# Adding extra points to include the different heights
# ----------------------------------------------------
# p20 == p2 
# p21 == p5 
# p22 == p3
# p23 == p4
# p24 == p8
# p25 == p9
# p26 == p14
# p27 == p15

pX.append(geompy.MakeVertex(fac*arr[2,0], fac*arr[2,1], fac*hghS))
pX.append(geompy.MakeVertex(fac*arr[5,0], fac*arr[5,1], fac*hghS))
pX.append(geompy.MakeVertex(fac*arr[3,0], fac*arr[3,1], fac*hghS))
pX.append(geompy.MakeVertex(fac*arr[4,0], fac*arr[4,1], fac*hghS))
pX.append(geompy.MakeVertex(fac*arr[8,0], fac*arr[8,1], fac*hghS))
pX.append(geompy.MakeVertex(fac*arr[9,0], fac*arr[9,1], fac*hghS))
pX.append(geompy.MakeVertex(fac*arr[14,0], fac*arr[14,1], fac*hghS))
pX.append(geompy.MakeVertex(fac*arr[15,0], fac*arr[15,1], fac*hghS))
                        
# ------------------------- #            
# Adding points to study
                      
point_ID = []
for i in range(int(len(pX)/2.)):
   point_ID.append(geompy.addToStudy(pX[i], "p"+repr(i)))
   gg.createAndDisplayGO(point_ID[i])
                            
name = "DropletFormation"

# ------------------------------------------------------------------------------- #
# Functions for creating different type of objects : BOX and TRAPEZOID

def CreateBoxTwoPoints(p1,p2):
    return geompy.MakeBoxTwoPnt(p1, p2)
    
def CreateTrapezoidEightPoints(p1, p2, p3, p4, p5, p6, p7, p8):
    tFace1 = geompy.MakeQuad4Vertices(pX[p1], pX[p2], pX[p3], pX[p4])
    tFace2 = geompy.MakeQuad4Vertices(pX[p5], pX[p6], pX[p7], pX[p8])
    return  geompy.MakeHexa2Faces(tFace1, tFace2)
            
def CreatePrismSixPoints(p1, p2, p3, p4, p5, p6):
    e1 = geompy.MakeEdge(pX[p1], pX[p2])
    e2 = geompy.MakeEdge(pX[p2], pX[p3])
    e3 = geompy.MakeEdge(pX[p3], pX[p1])
    e4 = geompy.MakeEdge(pX[p4], pX[p5])   
    e5 = geompy.MakeEdge(pX[p5], pX[p6])      
    e6 = geompy.MakeEdge(pX[p6], pX[p4])
    w1 = geompy.MakeWire([e1, e2, e3])
    w2 = geompy.MakeWire([e4, e5, e6])
    pFace1 = geompy.MakeQuad4Vertices(pX[p1], pX[p4], pX[p5], pX[p2])
    pFace2 = geompy.MakeQuad4Vertices(pX[p3], pX[p6], pX[p5], pX[p2])
    pFace3 = geompy.MakeQuad4Vertices(pX[p1], pX[p4], pX[p6], pX[p3])
    pFace4 = geompy.MakeFaces([w1], 0)
    pFace5 = geompy.MakeFaces([w2], 0)
    shell = geompy.MakeShell([pFace1, pFace2, pFace3, pFace4, pFace5])
    return geompy.MakeSolid([shell])

                                                            
# --------------------------------------------------------------------------------- #
# Creating boxes :
box1 = CreateBoxTwoPoints(pX[3], pX[nop + 1])
box2 = CreateBoxTwoPoints(pX[7], pX[nop + 5])
box3 = CreateBoxTwoPoints(pX[4], pX[2*nop])
box4 = CreateBoxTwoPoints(pX[9], pX[2*nop + 2])
box5 = CreateBoxTwoPoints(pX[2*nop], pX[nop + 4])  

# Creating Trapezoids
TSolid1 = CreateTrapezoidEightPoints(14, 15, 9, 8, 2*nop + 6, 2*nop + 7, 2*nop + 5, 2*nop + 4)
TSolid2 = CreateTrapezoidEightPoints(14, 16, 10, 8, 14 + nop, 16 + nop, 10 + nop, 8 + nop)
TSolid3 = CreateTrapezoidEightPoints(15, 17, 11, 9, 15 + nop, 17 + nop, 11 + nop, 9 + nop)
TSolid4 = CreateTrapezoidEightPoints(18, 16, 10, 12, 18 + nop, 16 + nop, 10 + nop, 12 + nop)
TSolid5 = CreateTrapezoidEightPoints(19, 17, 11, 13, 19 + nop, 17 + nop, 11 + nop, 13 + nop)
TSolid6 = CreateTrapezoidEightPoints(2*nop + 6, 2*nop + 7, 2*nop + 5, 2*nop + 4, nop + 14, nop + 15, nop + 9, nop + 8)

# Creating combined geometry

tmp = geompy.MakeFuse(box1, box2)
tmp = geompy.MakeFuse(tmp, box3)
tmp = geompy.MakeFuse(tmp, box4)
tmp = geompy.MakeFuse(tmp, box5)
tmp = geompy.MakeFuse(tmp, TSolid1)  
tmp = geompy.MakeFuse(tmp, TSolid2)
tmp = geompy.MakeFuse(tmp, TSolid3)  
tmp = geompy.MakeFuse(tmp, TSolid4)
tmp = geompy.MakeFuse(tmp, TSolid5)
final = geompy.MakeFuse(tmp, TSolid6)

# Adding Geometry to study

id_p1 = geompy.addToStudy(box1, "box1")
id_p2 = geompy.addToStudy(box2, "box2")
id_p3 = geompy.addToStudy(box3, "box3")
id_p4 = geompy.addToStudy(box4, "box4")
id_p5 = geompy.addToStudy(box5, "box5")
id_p6 = geompy.addToStudy(TSolid1, "TSolid1")
id_p7 = geompy.addToStudy(TSolid2, "TSolid2")
id_p8 = geompy.addToStudy(TSolid3, "TSolid3")
id_p9 = geompy.addToStudy(TSolid4, "TSolid4")
id_p10 = geompy.addToStudy(TSolid5, "TSolid5")
id_p11 = geompy.addToStudy(TSolid6, "TSolid6")
id_final = geompy.addToStudy(final, "Final") 

# Displaying Geometry
gg.createAndDisplayGO(id_final)

# Creating Groups : 2 Inlets, 2 Outlets and walls

# Adding Objects to groups

Inlet1 = geompy.CreateGroup(final, geompy.ShapeType["FACE"])
Inlet2 = geompy.CreateGroup(final, geompy.ShapeType["FACE"])
Outlet1 = geompy.CreateGroup(final, geompy.ShapeType["FACE"])
Outlet2 = geompy.CreateGroup(final, geompy.ShapeType["FACE"])
Walls = geompy.CreateGroup(final, geompy.ShapeType["FACE"])

# Adding Objects to respective groups 
SubFaceList = geompy.SubShapeAll(final, geompy.ShapeType["FACE"])

FaceID1 = geompy.GetSubShapeID(final, SubFaceList[0])
geompy.AddObject(Inlet1, FaceID1)     

FaceID2 = geompy.GetSubShapeID(final, SubFaceList[30])
geompy.AddObject(Inlet2, FaceID2)

FaceID3 = geompy.GetSubShapeID(final, SubFaceList[20])
geompy.AddObject(Outlet1, FaceID3)

FaceID4 = geompy.GetSubShapeID(final, SubFaceList[40])
geompy.AddObject(Outlet2, FaceID4)

# Adding all the objects to to Group "Walls" and removing the inlets/outlets

for i in range(len(SubFaceList)):
   FaceID = geompy.GetSubShapeID(final, SubFaceList[i])
   geompy.AddObject(Walls, FaceID)

geompy.RemoveObject(Walls, FaceID1)
geompy.RemoveObject(Walls, FaceID2)
geompy.RemoveObject(Walls, FaceID3)   
geompy.RemoveObject(Walls, FaceID4)

id_group1 = geompy.addToStudy(Inlet1, "Inlet1")
id_group2 = geompy.addToStudy(Inlet2, "Inlet2")
id_group3 = geompy.addToStudy(Outlet1, "Outlet1")
id_group4 = geompy.addToStudy(Outlet2, "Outlet2")
id_group5 = geompy.addToStudy(Walls, "Walls")
                                           
# --------------------------------------------------------------------------------- #
# MESHING

# Generating Tetrahedral Mesh
import smesh

# Creating Mesh
tetra = smesh.Mesh(final, name)

# Define 1D hypothesis      
algo1d = tetra.Segment()
algo1d.LocalLength(fac*0.02)

# Define 2D hypothesis
algo2d = tetra.Triangle()
algo2d.LengthFromEdges()

# Define 3D hypothesis
algo3d = tetra.Tetrahedron(smesh.NETGEN)
algo3d.MaxElementVolume(fac*0.006)

# Compute the mesh
tetra.Compute()

# Getting the information about the meshes : 

print "Information about mesh:"
print "Number of nodes       : ", tetra.NbNodes()
print "Number of edges       : ", tetra.NbEdges()
print "Number of faces       : ", tetra.NbFaces()
print "          triangles   : ", tetra.NbTriangles()
print "          quadrangles : ", tetra.NbQuadrangles()
print "          polygons    : ", tetra.NbPolygons()
print "Number of volumes     : ", tetra.NbVolumes()
print "          tetrahedrons: ", tetra.NbTetras()
print "          hexahedrons : ", tetra.NbHexas()
print "          prisms      : ", tetra.NbPrisms()
print "          pyramids    : ", tetra.NbPyramids()
print "          polyhedrons : ", tetra.NbPolyhedrons()

# Creating Mesh Groups from Geometry

sMesh1 = tetra.GroupOnGeom(Inlet1)
sMesh2 = tetra.GroupOnGeom(Inlet2)
sMesh3 = tetra.GroupOnGeom(Outlet1)
sMesh4 = tetra.GroupOnGeom(Outlet2)
sMesh5 = tetra.GroupOnGeom(Walls)

salome.sg.updateObjBrowser(1)

# Exporting to UNV File

tetra.ExportUNV("---- add path ---- /meshUNV.unv", 0)

# -----------------------------------------------end-of-file----------------------- #
