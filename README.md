# ThRend
Infrared rendering for thermography simulation

### Introduction

ThRend is a ray-tracing-based renderer of infrared radiation. Based on the Embree ray-tracing kernels, ThRend is designed to fastly generate simulated thermograms based on few input data. 
This software can be considered as a post-processing tool that takes the output of other thermal simulation software and allows to simulate the behavior of long-wave radiation reaching an infrared sensor (thermal camera).

### Input data

ThRend input data is handled through two configuration files: *viewSettings* and *materials*. 
*viewSettings* contains the configuration of the scene, camera and output images. 
*materials* describes the infrared properties of the materials to be used.












