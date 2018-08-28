#!/bin/sh

LIBS="-lann -lboost_system -lboost_filesystem -lboost_program_options -I/usr/local/include/opencv -I/usr/local/include -L/usr/local/lib -lopencv_ml -lopencv_objdetect -lopencv_shape -lopencv_stitching -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_imgproc -lopencv_flann -lopencv_core"

g++ -g  -o pers --std=c++11  *.cpp Wrappers/*.cpp SimPers/*.cpp Graphs/*.cpp GIComplex/*.cpp  $LIBS 
