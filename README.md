# ThRend
Infrared rendering for thermography simulation

### Introduction

ThRend is a ray-tracing-based renderer of infrared radiation. Based on the Embree ray-tracing kernels, ThRend is designed to fastly generate simulated thermograms based on few input data. 
This software can be considered as a post-processing tool that takes the output of other thermal simulation software and allows to simulate the behavior of long-wave radiation reaching an infrared sensor (thermal camera).

### Input data

ThRend input data is handled through two configuration files: *viewSettings* and *materials*. 
*viewSettings* contains the configuration of the scene, camera and output images. 
*materials* describes the infrared properties of the materials to be used. 

#### Example of *viewSettings* file: 
```
#scene settings                                               
sceneFile bayonne14hs.inp
skyTempsFile tsky  

# camera and image settings
# location of the camera:
cameraCenter  -87.6051 9.4538 1.0

# direction of the camera:
cameraDirection 0.9623 -0.21 0.172
cameraUp 0 0 1

# field of view of camera in vertical direction (in degrees)
fovVertical 24

# resolution of the image 
imageWidth 180
imageHeight 250

# Primary rays per pixel (for antialiasing)
aa 16

# Number of reflected rays per pixel
reflSamples 100;

# Colormap settings
colormapFile colormap
tmin 10
tmax 40

tmin_reflected -10
tmax_reflected 30
```
In this example, the scene geometry and nodal temperatures are loaded from the AVS UCD file *bayonne14hs.inp*, as indicated by the tag *sceneFile*. This kind of file can be exported from most thermal software, and its specifications can be found in [1]. For example, in Cast3m, you can export your CHPOINT *chp1* to this format with the following command: 
``` 
SORT AVS geo1 chp1
``` 
Currently, ThRend supports geometries (*geo1*) composed of quad and tri sufraces, but it should be very easy to adapt the code to handle other element types.

The *skyTempsFile* tag specifies the file where to find the sky temperatures. This file has only one line with 10 values indicating the temperature of the sky in different zenith angles (specified in kelvin, from the zenith into the horizon with 10 degrees steps):
``` tsky file
233.6 235.4 238.4 242.6 248 254.6 262.4 271.4 281.6 293.0
``` 
These temperatures are used for rays that do not hit any geometry. If you do not want to use this alternative, just input a fully closed mesh (eg. a room or a sky box).

The following tags set the camera properties, which are pretty self-explanatory. 
The tag *aa* indicates the number of primary rays per pixel for antialiasing. Please beware that the execution time of ThRend is linear with respect to this number, so try to keep it as small as possible. If you do not want antialiasing, just put *aa 1*.

The tag *reflSamples* indicates the number of reflected rays to cast per pixel. You want to have a number of rays that ensures sufficiently good sampling. 

The tag *colormapFile* indicates the file where to find the colormap specification of the output, using the same format as MATLAB colormaps (each line contains one RGB color). The tags *tmin* and *tmax* indicate the colormap temperature limits (in celsius) used for the output (for example, minimum temperature 10C and maximum 40C). *tmin_reflected* and *tmax_reflected* are the equivalent for the reflected temperature output image, which in this example are lower because most of the reflected temperature correspond to the cold sky. 

#### Example of *viewSettings* file: 
``` 
# material definition
name wood
UCD_id 2
normal_emissivity 0.95
diffuse_fraction 0.9
specular_lobe_size 200

name mortar
UCD_id 3
normal_emissivity 0.91
diffuse_fraction 0.85
specular_lobe_size 300

name glass
UCD_id 9
normal_emissivity 0.92
diffuse_fraction 0
specular_lobe_size 1.0E+5
``` 

This file contains the definition of the infrared properties of the materials. Each material is associated with a different color in the UCD file. The idea here is to have the elements grouped by color (e.g. Cast3m colors), where each color has its own material properties. The tag *name* indicates the name of the created material, while the tag *UCD_id* indicates the associated id in the geometry file.

The tag *normal_emissivity* indicates the emissivity value at normal direction. This value is used to generate a directional emissivity curve following Schlick approximation

### Execution

### Output data

### Source code and compilation

### References
[1] AVS-UCD file format description. site: https://dav.lbl.gov/archive/NERSC/Software/express/help6.1/help/reference/dvmac/UCD_Form.htm. Accessed: May 8, 2020.




