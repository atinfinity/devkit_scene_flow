# devkit_scene_flow
The development kit of Stereo Evaluation 2015(OpenCV version)

## Introduction
The original [development kit](http://kitti.is.tue.mpg.de/kitti/devkit_scene_flow.zip) of [Stereo Evaluation 2015](http://www.cvlibs.net/datasets/kitti/eval_scene_flow.php?benchmark=stereo) depends on [png++](http://savannah.nongnu.org/projects/pngpp/).
As a result, it is difficult to execute on Windows. So, I changed this development kit to execute on Windows and Linux.

## Requirements
- OpenCV 2.4.x or 3.x
- CMake

## Build Instructions
```
$ git clone https://github.com/atinfinity/devkit_scene_flow.git
$ cd devkit_scene_flow/devkit/cpp
$ mkdir build
$ cd build
$ cmake -D CMAKE_BUILD_TYPE=RELEASE ..
$ make
```
