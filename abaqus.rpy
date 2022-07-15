# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-5 replay file
# Internal Version: 2015_08_18-20.07.49 135153
# Run by Syngman on Fri Jul 15 09:39:12 2022
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=207.83332824707, 
    height=132.701400756836)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
execfile('create_geometry.py', __main__.__dict__)
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
p1 = mdb.models['Model-1'].parts['SPA']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
