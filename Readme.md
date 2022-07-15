# Multi-Chamber-Bending-SPA-Generator

This script is created for our final year project <strong>Soft Grippers For Automated Strawberry Harvesting</strong> soft gripper simulation and this is a modified version of research script of <a href="epfl.ch">EPFL</a> University

Multichamber bending soft pneumatic actuator generator for abaqus

How to use the script in windows OS

1.  Install Abaqus software and add it to windows path

2.  Install python 2.7 and add it to windows path

3.  Change script parameters according to your needs

    Select actuator = 'bnd' for bending soft pneumatic actuators and 'lin' for linear soft pneumatic actuators<br>

    Script is designed for ecoflex30. If you need to use an different material change the ogden3.mat file<br>

    Geometry parameters<br>

    Extrude a solid rectangle for body (assuming mirror-symmetry in x).<br>
    O = origin<br>
    H = body_H<br>
    W = body_W<br>
    L = body_L<br>

                  +-----W-----+
                 /           / |
               L           L   H
             /           /     |
           +-----W-----+       +
           |           |     /     Y
           |           H   L       ^    -Z
           INLET       | /         |   /
           O- - - - - -+           | /
               SYMM                O-----> X

    Extrude a cut for the inlet tunnel.<br>

    O = origin<br>
    H = body_H<br>
    h = inlet_H<br>
    L = body_L<br>
    l = wall<br>

             /                      /
           +----------L-----------+
           |                      |
           |                      | /
           |                      |
           H            +----l--- H
           |            h   CUT   |
           |            |         |
           +- - - - - - - - - - - O

4.  Open VSCODE in repo folder path or Open Command promt in repo folder path
5.  Run command <strong>abaqus cae script=create_geometry.py
