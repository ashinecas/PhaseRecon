function R = pmhRotationMatrix(phi,theta,psi)

%%**************************************************************************
%% System name:      Calibration from BB image
%% Module name:      pmhRotationMatrix.m
%% Version number:   1.0
%% Revision number:  00
%% Revision date:    3-Nov-2003
%%
%% 2003 (C) Copyright by Douglas J. Moseley
%%          Princess Margaret Hospital
%%
%%  Inputs:
%%
%%  Outputs:
%%
%%  Description:
%%      Compute rotation matrix
%%
%%      [ cos(theta)*cos(psi)  cos(phi)*cos(psi)
%%  Notes:
%%
%%*************************************************************************
%% References: 
%%
%%*************************************************************************
%% Revision History
%%	0.100    2003 11 3     Initial version
%%*************************************************************************

Theta = theta*pi/180;
Phi = phi*pi/180;
Psi = psi*pi/180;

R = [ 1           0         0
      0           cos(Phi)  sin(Phi)
      0          -sin(Phi)  cos(Phi)]  ...
    * ...
    [ cos(Theta)  0        -sin(Theta)
      0           1         0
      sin(Theta)  0         cos(Theta)] ...
    * ...
    [ cos(Psi)    sin(Psi)  0
     -sin(Psi)    cos(Psi)  0
      0           0         1         ];


return
   
