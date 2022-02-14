function [GC, MC, k1, k2, ABS, RMS] = implicitcurvature(TPMS,points)
%IMPLICTCURVATURE Calcualtes curvature based on the implicit TPMS functions
% Inputs:
%   TPMS - the TPMS implicit to be analysed
%   points - the list of points (verticies or face centroids) to evalue curvature at 
% Outputs:
%   GC - gaussian curvature (from determinant of the shape matrix)
%   MC - mean curvature (from divergance of unit normal)
%   k1 - maximum principle curvature (from GC and MC)
%   k2 - minimum principle curvature (from GC and MC)
%   ABS - absoloute curvature (from k1 and k2)
%   RMS - rms curvature (from GC and MC)
% Designed for use in the TPMS Designer Toolbox
%
%   Credit: Alistair Jones, 30/09/2021, RMIT University
%% Find Various Properties - Gaussian, Mean, Inclination

syms x y z;
assume(x,'real'); assume(y,'real'); assume(z,'real');

%Apply the scaling of the unit cell before differentiation
xyz = 2*pi*[x y z]./[TPMS.cellSize];

if TPMS.type=="network"
    U = TPMS.u;
else
    U = matlabFunction((TPMS.u(x,y,z)).^2);
end

% Create implicit functions for MC and GC using differential geometry
G = gradient(U(xyz(1),xyz(2),xyz(3)));
G2 = [gradient(G(1),[x y z]) gradient(G(2),[x y z]) gradient(G(1),[x y z])];
K = det([G2' G; G' 0]);
H = -1/2*divergence(G/norm(G));
MCf = matlabFunction(H,'Vars',[x y z]);
GCf = matlabFunction(-K/(G'*G).^2,'Vars',[x y z]);


%Rotations and phase shifts are shape preserving - these can be applied by
%inversely transforming the points of interest.
rZ = [cosd(TPMS.Rxyz(3)) -sind(TPMS.Rxyz(3)) 0; sind(TPMS.Rxyz(3)) cosd(TPMS.Rxyz(3)) 0; 0 0 1];
rY = [cosd(TPMS.Rxyz(2)) 0 -sind(TPMS.Rxyz(2)); 0 1 0; sind(TPMS.Rxyz(2)) 0 cosd(TPMS.Rxyz(2))];
rX = [1 0 0; 0 cosd(TPMS.Rxyz(1)) -sind(TPMS.Rxyz(1)); 0 sind(TPMS.Rxyz(1)) cosd(TPMS.Rxyz(1))];
trans = [TPMS.offset(1) TPMS.offset(2) TPMS.offset(3)];
pad = [0; 0; 0]; T = rX*rY*rZ;
tform=affine3d([T pad; trans 1]);
points = transformPointsInverse(tform,points);

%Calculate the properties at the transformed points of interest
GC = GCf(points(:,1),points(:,2),points(:,3));
MC = MCf(points(:,1),points(:,2),points(:,3));
k1 = real(MC+sqrt(MC.^2-GC));
k2 = real(MC-sqrt(MC.^2-GC));
ABS = abs(k1)+abs(k2);
RMS = real(sqrt(4*MC.^2-2*GC));
end

