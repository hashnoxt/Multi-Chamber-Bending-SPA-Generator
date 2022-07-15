# -*- coding: mbcs -*-
# -*- coding: utf-8 -*-

import csv
import re
from math import *
from pyGene4SPA.gene import FloatGene
from pyGene4SPA.organism import Organism
from pyGene4SPA.population import Population
from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import *
from driverUtils import executeOnCaeStartup
import model_utility as MU
import sys
import pprint
import os
import os.path
import shutil
sys.path.append('/usr/lib/python2.7/')
import argparse
import inspect
import numpy as np
import math

# return the model given a string model name


def get_abq_model(mname):
    mn = mname.lower()
    if mn == 'ab':
        return ARRUDA_BOYCE
    if mn == 'mr':
        return MOONEY_RIVLIN
    if mn == 'neoh':
        return NEO_HOOKE
    if mn == 'ogden':
        return OGDEN
    if mn == 'poly':
        return POLYNOMIAL
    if mn == 'rpoly':
        return REDUCED_POLYNOMIAL
    if mn == 'vdw':
        return VAN_DER_WAALS
    if mn == 'yeoh':
        return YEOH
    print 'ERROR: invalid model name \'' + mname + '\''
    exit(1)
# return the max result from the field output

def max_result(result,o):
    result_field, result_invariant = result
    __max = 0
    frames = o.steps['Step-1'].frames
    f1 = frames[-1]  # Acquire the last frame
    fop = f1.fieldOutputs
    if(fop.has_key(result_field)):
        measureSet = fop[result_field]
        for measureValue in measureSet.values:
            if result_invariant:
                if hasattr(measureValue, result_invariant.lower()):
                    val = getattr(measureValue, result_invariant.lower())
                else:
                    raise ValueError(
                        'Field value does not have invariant %s' % (result_invariant,))
            else:
                val = measureValue.data
            if(val > __max):
                __max = val
    else:
        raise ValueError('Field output does not have field %s' %
                         (results_field,))
    return __max


mat = MU.read_matfile('ogden3.mat')




def performCAE(polist):


    quadelms = True
    ALE = False
    time = 200

    nc = 10  # number of chambers
    wall = 8
    cwall = 3*wall/4
    gap = wall/3
    uplayer = 4
    btlayer = 8
    #chamber_size
    chamber_L = polist[0]
    chamber_H = polist[1]
    chamber_W = chamber_L

    offset = 1.0  # extended for length
    Inexlayer = 0.35*btlayer
    inlet_H = (btlayer - Inexlayer) / 2
    inlet_W = chamber_W / 3

    body_H = chamber_H + uplayer + btlayer
    body_W = chamber_W + (wall/2)
    body_L = 2 * (wall + offset) + nc * chamber_L + (nc - 1) * cwall

    mass_scaling = 1.0
    film_thickness = 1.0
    mesh_size = 4.0
    #pressure 
    PRESSURE_C = 0.03
    Beta = 0.95
    PRESSURE_A = (math.log(PRESSURE_C+1, 5))/30
    PRESSURE_B = PRESSURE_A*Beta + PRESSURE_C*(1-Beta)


    # Open the viewport
    #myViewport = session.Viewport(
    #    name='HXZSPA', origin=(0, 0), width=250, height=150)
    executeOnCaeStartup()
    Mdb()

    # create the part
    s = mdb.models['Model-1'].ConstrainedSketch(name='profile', sheetSize=200.0)
    s.setPrimaryObject(option=STANDALONE)

    # Extrude a solid rectangle for body (assuming mirror-symmetry in x).
    s.rectangle(point1=(0.0, 0.0), point2=(body_W, body_H))
    p = mdb.models['Model-1'].Part(name='SPA',
                                   dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['SPA']
    p.BaseSolidExtrude(sketch=s, depth=body_L)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['profile']

    # Extrude cuts for the air inlet tunnel
    f = p.faces
    e = p.edges
    face = f.findAt(((0.0, 0.1, 0.1),),)[0]
    edge = e.findAt(((0.0, body_H - 0.1, body_L),),)[0]
    t = p.MakeSketchTransform(sketchPlane=face, sketchUpEdge=edge,
                              sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, body_L))
    s1 = mdb.models['Model-1'].ConstrainedSketch(
        name='__profile__', sheetSize=45.6, gridSpacing=1.14, transform=t)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    s1.rectangle(point1=(0.0, btlayer - inlet_H),
                 point2=(-body_L + wall + offset, btlayer))
    p.CutExtrude(sketchPlane=face, sketchUpEdge=edge, sketchPlaneSide=SIDE1,
                 sketchOrientation=RIGHT, sketch=s1, depth=inlet_W, flipExtrudeDirection=OFF)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']

    # Extrude cuts for the air chambers
    f = p.faces
    e = p.edges
    face = f.findAt(((0.0, 0.1, 0.1),),)[0]
    edge = e.findAt(((0.0, body_H - 0.1, body_L),),)[0]
    t = p.MakeSketchTransform(sketchPlane=face, sketchUpEdge=edge,
                              sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, body_L))
    s = mdb.models['Model-1'].ConstrainedSketch(
        name='profile', sheetSize=45.6, gridSpacing=1.14, transform=t)
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    for i in range(nc):
      xoff = offset + wall + i * (cwall + chamber_L)
      s.rectangle(point1=(-xoff, btlayer),
                  point2=(-xoff - chamber_L, btlayer + chamber_H))

    p.CutExtrude(sketchPlane=face, sketchUpEdge=edge, sketchPlaneSide=SIDE1,
                 sketchOrientation=RIGHT, sketch=s, depth=chamber_W, flipExtrudeDirection=OFF)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['profile']

    # Extrude cuts for the gaps between chambers
    f = p.faces
    e = p.edges
    face = f.findAt(((0.0, 0.1, 0.1),),)[0]
    edge = e.findAt(((0.0, body_H - 0.1, body_L),),)[0]
    t = p.MakeSketchTransform(sketchPlane=face, sketchUpEdge=edge,
                              sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, body_L))
    s = mdb.models['Model-1'].ConstrainedSketch(
        name='profile', sheetSize=45.6, gridSpacing=1.14, transform=t)
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    for i in range(nc - 1):
      xoff = offset + wall + chamber_L + i * (cwall + chamber_L)
      s.rectangle(point1=(-xoff-((cwall-gap)/2), body_H),
                  point2=(-xoff -((cwall-gap)/2)-gap, body_H - chamber_H*1.2))

    p.CutExtrude(sketchPlane=face, sketchUpEdge=edge, sketchPlaneSide=SIDE1,
                 sketchOrientation=RIGHT, sketch=s, depth=body_W, flipExtrudeDirection=OFF)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['profile']

    # Define some variables
    d = p.datums
    c = p.cells
    f = p.faces
    e = p.edges
    v = p.vertices

    # Create sets and surfaces
    p.Set(cells=c, name='ALL')
    # Symmetry boundary condition in width.
    faces = f.findAt(((0.0, body_H-0.1, 0.1),),)
    p.Set(faces=faces, name='SYMM_W')
    # Pressure set, for measuring purposes
    zoff = body_L - offset - wall - (nc - 1) * (cwall + chamber_L)
    # chamber-left
    faces = f.findAt(
        ((chamber_W - 0.1, btlayer + chamber_H - 0.1, zoff - chamber_L),),)
    p.Set(faces=faces, name='PRESSURE')
    # Surface for measuring displacements (same for linear or bending).
    faces = f.findAt(((0.1, 0.1, 0.0),),)   # back side.
    p.Surface(side1Faces=faces, name='MEASURED')

    # Partition the part into two regions
    pickedCells = c.findAt(((0.1, 0.1, 0.1),))
    p.DatumPlaneByPrincipalPlane(
        principalPlane=XZPLANE, offset=Inexlayer)
    p.PartitionCellByDatumPlane(datumPlane=d.values()[-1], cells=pickedCells)
    # Split the middle section along the back of the chambers.
    pickedCells = c.findAt(((0.1, Inexlayer+0.1, 0.1), ))
    #pickedCells = c.getSequenceFromMask(mask=('[#2 ]', ), )
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=chamber_W)
    p.PartitionCellByDatumPlane(datumPlane=d.values()[-1], cells=pickedCells)
    # Partition the membrane for the tube BC, if necessary.
    # pickedFaces = f.findAt(((0.1, 0.0, wall + offset + 0.1), ))  # Membrane.
    # p.DatumPlaneByPrincipalPlane(
    # principalPlane=XYPLANE, offset=body_L - wall - offset)
    # p.PartitionFaceByDatumPlane(datumPlane=d.values()[-1], faces=pickedFaces)

    # Additional sets and surfaces after partitioning

    # Tube(fixed) boundary condition
    faces = f.findAt(((body_W - 0.1, 0.1 + btlayer, body_L),),  # extrenal face 1
                     ((inlet_W + 0.1, Inexlayer + 0.1, body_L),),  # external face 2
                     ((0.1, 0.1, body_L),),)  # external face 3
    p.Set(faces=faces, name='TUBE')

    # pressure boundary condition in chambers and connecting tubes
    #Pressure A
    for i in range(nc):
      zoff = wall + offset + i * (cwall + chamber_L)
      faces = f.findAt(((chamber_W - 0.1, btlayer + chamber_H, zoff + 0.1),),  # chamber top
                       ((chamber_W - 0.1, btlayer, zoff + 0.1),),  # chamber bottom
                       ((chamber_W - 0.1, btlayer + 0.1, zoff),),  # chamber left
                       ((chamber_W - 0.1, btlayer + 0.1, zoff + chamber_L),), # chamber right
                       ((chamber_W, btlayer + 0.1, zoff + 0.1),),)  # chamber down

      faces2 = f.findAt(
          ((0.1, btlayer, zoff + chamber_L + 0.1),),)  # inlet-top
      if i == 0:
        allfaces = faces + faces2
      else:
        allfaces = allfaces + faces + faces2

    p.Surface(side1Faces=allfaces, name='PRESSURE_A')

    # For a bending actuator, also apply the pressure to the air inlet tube
    #allfaces = p.surfaces['PRESSURE'].faces
    #allfaces = allfaces + \
    #    f.findAt(((inlet_W, btlayer - 0.1, wall + offset + 0.1),))
    #allfaces = allfaces + \
    #    f.findAt(((0.1, btlayer - inlet_H, wall + offset + 0.1),))
    #p.Surface(side1Faces=allfaces, name='PRESSURE')

    #Pressure B
    facesB = f.findAt(((chamber_W-0.1,btlayer+0.1,wall+offset),))
    facesB = facesB + f.findAt(((inlet_W, btlayer - 0.1, wall + offset + 0.1),))
    p.Surface(side1Faces=facesB,name='PRESSURE_B')

    #pressure C
    facesC = f.findAt(((0.1, btlayer - inlet_H, wall + offset + 0.1),))
    p.Surface(side1Faces=facesC,name='PRESSURE_C')

    # Create the material for the Inextensible layer
    mdb.models['Model-1'].Material(name='SILK')
    mdb.models['Model-1'].materials['SILK'].Density(
        table=((mass_scaling * 0.0013, ), ))
    mdb.models['Model-1'].materials['SILK'].Elastic(table=((4100.0, 0.49), ))
    # Create a material section. Note that increasing the thickness will increase
    # the tension resistance, but will take longer to converge.
    # mdb.models['Model-1'].MembraneSection(name='SILK', material='SILK', thicknessType=UNIFORM,
                                          # thickness=film_thickness, thicknessField='', poissonDefinition=DEFAULT)
    mdb.models['Model-1'].HomogeneousSolidSection(
        name='SILK', material='SILK', thickness=film_thickness)

    cells = c.findAt(((0.1,0.1,0.1),))
    region1 = p.Set(cells=cells, name='SILK_LAYER')
    p.SectionAssignment(region=region1, sectionName='SILK', offset=0.0,
        offsetType=MIDDLE_SURFACE, offsetField='',
        thicknessAssignment=FROM_SECTION)


    # Set the output frequency
    output_frequency = 20

    #Create the material for the Extensible layer
    if mat['density']==0.0: print 'WARNING: material density equals 0.0'
    if mat['matname']=='???': print 'WARNING: material name is given as',mat['matname']
    mdb.models['Model-1'].Material(name=mat['matname'])
    mdb.models['Model-1'].materials[mat['matname']].Density(table=((mass_scaling*mat['density'], ), ))
    allparams = mat['params']+mat['D']
    if 'order' in mat:
        mdb.models['Model-1'].materials[mat['matname']].Hyperelastic(
            materialType=ISOTROPIC, testData=OFF, volumetricResponse=VOLUMETRIC_DATA,
            type=get_abq_model(mat['model']), n=mat['order'], table=(allparams, ),
            moduliTimeScale=LONG_TERM)
            # moduliTimeScale=INSTANTANEOUS)
    else:
        mdb.models['Model-1'].materials[mat['matname']].Hyperelastic(
            materialType=ISOTROPIC, testData=OFF, volumetricResponse=VOLUMETRIC_DATA,
            type=get_abq_model(mat['model']), table=(allparams, ),
            moduliTimeScale=LONG_TERM)
            # moduliTimeScale=INSTANTANEOUS)
    # Support hysteretic properties.
    if 'hyst_S' in mat:
        mdb.models['Model-1'].materials[mat['matname']].hyperelastic.Hysteresis(
            table=((mat['hyst_S'], mat['hyst_A'], mat['hyst_m'], mat['hyst_C']), ))
    mdb.models['Model-1'].HomogeneousSolidSection(name='SPA',material=mat['matname'],thickness=None)

    #Collect cell for the extensible layer
    cella = c.findAt(((0.1,Inexlayer+0.1,0.1),))
    cellb = c.findAt(((body_W-0.1,Inexlayer+0.1,0.1),))
    allcells = cella+cellb
    region2 = p.Set(cells=allcells,name='SPA_LAYER')
    p.SectionAssignment(region=region2, sectionName='SPA', offset=0.0,
        offsetType=MIDDLE_SURFACE, offsetField='',
        thicknessAssignment=FROM_SECTION)

    #Mesh the geometry

    #Using Tetrahedron and hexahedron
    #p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    #pickedRegions = c.findAt(((0.1,body_H-0.1,0.1),))
    #p.setMeshControls(regions=pickedRegions,elemShape=TET,technique=FREE)
    #elemType1 = mesh.ElemType(elemCode=C3D20RH, elemLibrary=STANDARD)
    #elemType2 = mesh.ElemType(elemCode=C3D15H, elemLibrary=STANDARD)
    #elemType3 = mesh.ElemType(elemCode=C3D10H, elemLibrary=STANDARD)
    #region = p.sets['ALL']
    #p.setElementType(regions=region, elemTypes=(elemType1, elemType2, elemType3))
    #p.generateMesh()

    #Only Tetrahedron
    pickedRegions = c.findAt(((0.1, body_H-0.1, 0.1), ), ((body_W-0.1, body_H-0.1, 
        0.1), ), ((0.1, 0.1, 0.1), ))
    p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
    cells = c.findAt(((0.1, body_H-0.1, 0.1), ), ((body_W-0.1, body_H-0.1, 
        0.1), ), ((0.1, 0.1, 0.1), ))
    pickedRegions =(cells, )
    elemType1 = mesh.ElemType(elemCode=C3D20RH, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D15H, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D10H, elemLibrary=STANDARD)
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2,elemType3))
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    # Create the assembly
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name='SPA',part=p,dependent=ON)

    #Create a step
    inc = min(15, 0.1*time)
    maxnuminc=100000
    #bnd actuator
    init_inc = inc/10.0
    mdb.models['Model-1'].ImplicitDynamicsStep(name='Step-1', previous='Initial',
                timePeriod=time, application=QUASI_STATIC, initialInc=init_inc, minInc=1e-08, maxInc=inc,
                maxNumInc=maxnuminc, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, initialConditions=OFF, nlgeom=ON)
    print 'No damping'

    if ALE:
      region = a.instances['SPA'].sets['ALL']
      mdb.models['Model-1'].steps['Step-1'].AdaptiveMeshDomain(region=region,controls=None,frequency=100,
        meshSweeps=25)   #Default=10,1
     
    #Create the loads and boundary conditions
    region = a.instances['SPA'].sets['SYMM_W']
    mdb.models['Model-1'].XsymmBC(name='SYMM_W',createStepName='Step-1',region=region,localCsys=None)
    region = a.instances['SPA'].sets['TUBE']
    mdb.models['Model-1'].DisplacementBC(name='TUBE', createStepName='Initial',
            region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    #pressure loads
    region = a.instances['SPA'].surfaces['PRESSURE_A']
    mdb.models['Model-1'].Pressure(name='PRESSURE_A', createStepName='Step-1',
                region=region, distributionType=UNIFORM, field='', magnitude=PRESSURE_A,
                amplitude=UNSET)

    region = a.instances['SPA'].surfaces['PRESSURE_B']
    mdb.models['Model-1'].Pressure(name='PRESSURE_B', createStepName='Step-1',
                region=region, distributionType=UNIFORM, field='', magnitude=PRESSURE_B,
                amplitude=UNSET)

    region = a.instances['SPA'].surfaces['PRESSURE_C']
    mdb.models['Model-1'].Pressure(name='PRESSURE_C', createStepName='Step-1',
                region=region, distributionType=UNIFORM, field='', magnitude=PRESSURE_C,
                amplitude=UNSET)

    # Create a reference point for calculating the length of the actuator.
    a.ReferencePoint(point=(body_W/2.0, body_H/2.0, body_L))
    rp = a.referencePoints[a.referencePoints.keys()[0]]
    a.Set(referencePoints=(rp,), name='LPT')
    region1 = a.sets['LPT']
    region2 = a.instances['SPA'].sets['TUBE']
    mdb.models['Model-1'].Coupling(name='LENGTH_MEASUREMENT', controlPoint=region1,
            surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING,
            weightingMethod=UNIFORM, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON,
            ur2=ON, ur3=ON)

    # Adjust the output requests.
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(frequency=output_frequency,variables=(
            'S',                          # Stress components and invariants.
            # 'PE', 'PEEQ', 'PEMAG',        # Plastic strain, Equivalent plastic strain, Plastic strain magnitude.
            'EE', 'IE', 'NE', 'LE',       # Elastic strain, Inelastic strain, Nominal strain, Logarithmic strain.
            'U', 'V', 'A',                # Displacement, Velocity, Acceleration.
            'RF', 'CF', 'P',              # Reaction forces and moments, Concentrated forces and moments, Pressure loads.
            'CSTRESS', 'CDISP', 'CFORCE', # Contact stresses, Contact displacements, Contact forces.
            'ENER',                       # All energy magnitudes.
            # 'EVOL',                       # Element volume.
            ))
    regionDef = a.allInstances['SPA'].sets['PRESSURE']
    mdb.models['Model-1'].FieldOutputRequest(name='PRESSURE', createStepName='Step-1',
            variables=('P', ), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, frequency=output_frequency)
    regionDef=a.sets['LPT']
    mdb.models['Model-1'].HistoryOutputRequest(name='LENGTH', variables=(('COOR3'),),
            createStepName='Step-1', region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, frequency=output_frequency)
        # H-Output-1 is the default history output request with the default variables.
    mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(frequency=output_frequency)

    # For a displacement test, add a nodeset to measure displacements.
    verts = v.findAt(((0.0, 0.0, 0.0), ))
    p.Set(vertices=verts, name='UPT')
    # Create a history output request for the displacement at the RP.
    regionDef = a.allInstances['SPA'].sets['UPT']
    mdb.models['Model-1'].HistoryOutputRequest(name='DISPLACEMENT', frequency=output_frequency,
                createStepName='Step-1', region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE,
                variables=('U1', 'U2', 'U3', 'COOR3'))
    #Create a job
    jname = 'MySPA_Analysis'
    myjob=mdb.Job(name=jname, model='Model-1', description='', type=ANALYSIS,
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
            explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF,
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
            scratch='', multiprocessingMode=DEFAULT, numCpus=8, numDomains=8,
            numGPUs=0)
    myjob.submit(consistencyChecking=OFF)
    myjob.waitForCompletion()
    #Extracts the maximum result from the Abaqus odb
    o = session.openOdb(name=jname+'.odb',readOnly=True)
    max_mises = max_result(('U','magnitude'),o)
    Result = max_mises
    

    return Result


def main():
    i = 0
    outFile = csv.writer(file('out_perform.csv', 'wb'))

    polist=[5.5,7.5]

    i=performCAE(polist)

    outStr=['Result',str(i)]
    outFile.writerow(outStr)


if __name__ == '__main__':
    main()
     