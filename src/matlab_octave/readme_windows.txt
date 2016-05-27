In order to use Optizelle in MATLAB, we must set the path using the following steps:

1.  In the MATLAB console, type

    userpath

2.  At the location found in step 1, add or modify the file "startup.m" with the line

    addpath('[INSTALL_ROOT]\share\optizelle\matlab')

where "[INSTALL_ROOT]" denotes the Optizelle install location.  Typically, this is something like

    C:\Program Files\Optizelle
