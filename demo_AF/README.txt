%%%%%%%%%%%%%   Run this under a Linux computer %%%%%%%%%%%%%%%

%% Unzip the test.zip file into the same directory as the C code.

%% Compile C code to executable and outputs with a.out 
gcc ExampleCcode.c -mcmodel=large -lm

%% Run the executable to produce *.vtk files for outputs with each time step
./a.out

%% Use free open source software "Paraview" to open the *.vtk outputs (https://www.paraview.org/download/) and the associated tutorials are also available on Paraview's website.


