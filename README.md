# ThRend
Infrared rendering for thermography simulation

### Introduction

ThRend is a ray-tracing-based renderer of infrared radiation. Based on the Embree ray-tracing kernels, ThRend is designed to fastly generate simulated thermograms based on few input data. 
This software can be considered as a post-processing tool that takes the output of other thermal simulation software and allows to simulate the behavior of long-wave radiation reaching an infrared sensor (thermal camera).

### Input data

ThRend input data is handled through two configuration files: *viewSettings* and *materials*. 
*viewSettings* contains the configuration of the scene, camera and output images. 
*materials* describes the infrared properties of the materials to be used. 

Example *viewSettings* file: 
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
In this example, the scene geometry and nodal temperatures are loaded from the AVS UCD file bayonne14hs.inp. This kind of file can be exported from most thermal software, and its specifications can be found in [1]. For example, in Cast3m, you can export your CHPOINT *chp1* to this format with the following command: 
``` 
SORT AVS geo1 chp1
``` 
Currently, ThRend supports geometries (*geo1*) composed of quad and tri sufraces. 




# References
[1] AVS-UCD file format description. site: https://dav.lbl.gov/archive/NERSC/Software/express/help6.1/help/reference/dvmac/UCD_Form.htm. Accessed: May 8, 2020.




