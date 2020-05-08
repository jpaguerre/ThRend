# ThRend
**Infrared rendering for thermography simulation**
<p align="center">
<img src="https://github.com/jpaguerre/ThRend/blob/master/README-IMGS/summary.png" width="90%" alt="centered image">
</p>

### Introduction

ThRend is a ray-tracing-based renderer of infrared radiation. Based on the Embree ray-tracing kernels, ThRend is designed to rapidly generate simulated thermograms based on few input data. 
This software can be considered as a post-processing tool that takes the output of other thermal simulation software and allows to simulate the behavior of long-wave radiation reaching an infrared sensor (thermal camera). ThRend gives the user the possibility of trying different emissivity and reflectivity configurations to render thermal images.

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
Currently, ThRend supports geometries (*geo1*) composed of quad and tri surfaces, but it should be very easy to adapt the code to handle other element types.

The *skyTempsFile* tag specifies the file where to find the sky temperatures. This file has only one line with 10 values indicating the temperature of the sky in different zenith angles (specified in kelvin, from the zenith into the horizon with 10 degrees steps):
``` tsky file
233.6 235.4 238.4 242.6 248 254.6 262.4 271.4 281.6 293.0
``` 
These temperatures are used for rays that do not hit any geometry. If you do not want to use this alternative, just input a fully closed mesh (eg. a room or a sky box).

The following tags set the camera properties, which are pretty self-explanatory. 
The tag *aa* indicates the number of primary rays per pixel for antialiasing. Please beware that the execution time of ThRend is linear with respect to this number, so try to keep it as small as possible. If you do not want antialiasing, just put *aa 1*.

The tag *reflSamples* indicates the number of reflected rays to cast per pixel. You want to have a number of rays that ensures sufficiently good sampling. 

The tag *colormapFile* indicates the file where to find the colormap specification of the output, using the same format as MATLAB colormaps (each line contains one RGB color). The tags *tmin* and *tmax* indicate the colormap temperature limits (in Celsius) used for the output (for example, minimum temperature 10C and maximum 40C). *tmin_reflected* and *tmax_reflected* are the equivalent for the reflected temperature output image, which in this example are lower because most of the reflected temperature correspond to the cold sky. 

#### Example of *materials* file: 
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

The tag *normal_emissivity* indicates the emissivity value at normal direction &epsilon;<sub>n</sub>. This value is used to generate a directional emissivity curve following Schlick's approximation [2]:

&epsilon; (&theta;) =  &epsilon;<sub>n</sub> - &epsilon;<sub>n</sub>(1 - cos(&theta;))<sup>5</sup>

The tag *diffuse_fraction* allows to generate an interpolation between the Schlick (Fresnel) emissivity curve and a diffuse curve of constant emissivity (using the value defined at *normal_emissivity*). This can be used to describe different levels of roughness in the materials, where for example glass is 0% diffuse, while wood is 90% diffuse. Note that this values define the emissivity behavior only; the reflectivity is handled with the next tag. Some examples of interpolated curves are:

![Emissivities](https://github.com/jpaguerre/ThRend/blob/master/README-IMGS/emissivity.png)

Two types of reflections can be used in ThRend: glossy and diffuse reflections. 
Glossy reflections are handled through importance sampling of the modified Phong BRDF model [3]. Hammersley sampling [4] is used to choose the ray directions. The tag *specular_lobe_size* indicates the size of the specular lobe; some examples of this value are:

<p align="center">
<img src="https://github.com/jpaguerre/ThRend/blob/master/README-IMGS/specsize.png" width="65%" alt="centered image">
</p>

See that a big value of *specular_lobe_size* implies a specular reflection behavior, this is why, in the example above, the material *glass* has a lobe size of 1.0E+5.

Diffuse reflections, on the other hand, are handled with Beckers-and-Beckers view factor sampling [5]. A purely diffuse reflection can be obtained by setting *specular_lobe_size* to be -1:
``` 
name diffuseMaterialName
UCD_id 11
normal_emissivity 0.9
diffuse_fraction 0.7
specular_lobe_size -1
``` 

One final note about material definition is that ThRend supports the definition of custom emissivity curves, such as in this example:
``` 
#custom material, defined with a curve of 91 emissivity values (from 90 to 0 degrees of viewing angle)
name customMat1
UCD_id 12
emissivity_curve 0.4650 0.5115 0.5528 0.5895 0.6223 0.6517 0.6779 0.7015 0.7228 0.7419 0.7591 0.7747 0.7889 0.8016 0.8132 0.8238 0.8333 0.8420 0.8500 0.8572 0.8637 0.8697 0.8752 0.8802 0.8847 0.8889 0.8927 0.8961 0.8993 0.9022 0.9048 0.9072 0.9095 0.9115 0.9133 0.9150 0.9165 0.9179 0.9192 0.9204 0.9215 0.9224 0.9233 0.9241 0.9249 0.9255 0.9261 0.9267 0.9272 0.9277 0.9281 0.9284 0.9288 0.9291 0.9294 0.9296 0.9298 0.9300 0.9302 0.9303 0.9305 0.9306 0.9307 0.9308 0.9309 0.9310 0.9310 0.9311 0.9312 0.9312 0.9312 0.9313 0.9313 0.9313 0.9313 0.9313 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 0.9314 
specular_lobe_size 200
``` 
The tag *emissivity_curve* tag replaces the *normal_emissivity* and *diffuse_fraction* tags. It allows defining a custom curve by entering 91 emissivity values: one value for each incidence angle (from 90 to 0 degrees).

### Execution
If you want to try ThRend without compiling it, just download the folder "executable/". You will find all the configuration files and an example of scene. Go into "binary/", and execute "ThRend.exe". You should see a console similar to this one:

<p align="center">
<img src="https://github.com/jpaguerre/ThRend/blob/master/README-IMGS/console.png" width="45%" alt="centered image">
</p>

The results will be saved in the folder "executable/results". This executable is only compatible with Windows 8/10 x64. If you need a Linux executable, do not hesitate to contact me.

### Output data
ThRend generates one main result and 4 auxiliary files. The main result is the file "apparent.png", which stores the rendered apparent surface temperature of the scene. The 4 auxiliary files are: "real.png", which shows the result as if every material was a blackbody, "emis.png", which stores a grayscale image showing the computed emissivity for each pixel, "refl.png" which shows the reflected temperature for each pixel (using the second colormap scale defined in *viewSettings*), and the file "temps", which stores the apparent temperatures as a matrix of numbers that can be loaded directly into Matlab. 

This is the output of the default project that comes with the executable:

<p align="center">
<img src="https://github.com/jpaguerre/ThRend/blob/master/README-IMGS/results.png" width="70%" alt="centered image">
</p>

### Source code and compilation
ThRend was developed in C++ with Visual Studio 2013. All the source code and VS projects are uploaded in this Github page. There are many dependencies but they are all portable, so you should be able to download and compile the project directly. The only prerequisite is Visual Studio 2013 or greater. If you need a Linux compilation project, do not hesitate to contact me. 

The libraries used by ThRend are the following:
1. Intel Embree 3 for ray tracing operations (https://www.embree.org/)
2. GLM for vector and matrix operations (https://glm.g-truc.net/0.9.9/index.html).
3. FreeImage for saving images (http://freeimage.sourceforge.net/).

### References
[1] AVS-UCD file format description. site: https://dav.lbl.gov/archive/NERSC/Software/express/help6.1/help/reference/dvmac/UCD_Form.htm. Accessed: May 8, 2020.

[2] Schlick, C. (1994, August). An inexpensive BRDF model for physically‐based rendering. In Computer graphics forum (Vol. 13, No. 3, pp. 233-246). Edinburgh, UK: Blackwell Science Ltd.

[3] Lafortune, Eric and Willems, Yves. Using the modified Phong reflectance model for physically based rendering. Report 197. Departement Computerwetenschappen, KU Leuven, Nov. 1994, 19 6.

[4] Suffern, K. (2016). Ray Tracing from the Ground up. CRC Press.

[5] Beckers, B., and Beckers, P. (2016, September). Fast and accurate view factor generation. In FICUP, An International Conference on Urban Physics (Vol. 9).



