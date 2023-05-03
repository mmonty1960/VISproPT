# VISproPT
VISproPT is the innovative instrument created at ENEA-Casaccia (Italy) for measuring the 3D shape of parabolic-trough (PT) panels for CSP solar collectors.

This C++ software covers all the steps of the 3D shape measurement of parabolic-trough reflective panels acchomplished with the VISproPT instrument:
    1) Camera-lens calibration, for image undistortion
    2) Instrument calibration
    3) Image-processing for evaluating: i) 3D shape (slopes dz/dx, dz/dy and height z), ii) deviations from the ideal shape and, last but not least, iii) evaluation of the intercept factor at a given longitudinal angle.
    
Installation on Linux:
  - install the libraries: Qt, OpenCV, Cminpack, cblas, blas, glib-2.0
  - download the source files (main.cpp, VISproPT.h, VISproPT.cpp, VISproPT.ui and VISproPT.pro) in a new folder in your PC
  - open a terminal window and change directory to the folder containig the source files
  - compile the project with the command "qmake"
  - complete with the command "make"
  
Process the exemplary images: 
  - Unzip the file RGinterior#58.7z contains the sequence of images acquired with the two cameras of VISproPT for an exemplary PT inner panel for LS3 solar collectors, as well as the configuration txt file containing the main instrument parameters
  - launch the software by a terminal window
  - load the configuration file
  - surf to the "Analysis" tab
  - push the button "Process!" and observe the progress in the terminal
  - at the process completion several countormaps will be drawn 
