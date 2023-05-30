function [GC, MC, k1, k2] = curvatureImplicit(TPMS,points)
%IMPLICTCURVATURE Calcualtes curvature based on the implicit TPMS functions
% Inputs:
%   TPMS - the TPMS implicit to be analysed
%   points - the list of points (verticies or face centroids) to evalue curvature at 
% Outputs:
%   GC - gaussian curvature (from determinant of the shape matrix)
%   MC - mean curvature (from divergance of unit normal)
%   k1 - maximum principle curvature (from GC and MC)
%   k2 - minimum principle curvature (from GC and MC)
% Designed for use in the TPMS Designer Toolbox
%
%   Credit: Alistair Jones, 30/09/2021, RMIT University
%% Find Various Properties - gaussian curvature, mean curvature, k1, k2

syms x y z;
assume(x,'real'); assume(y,'real'); assume(z,'real');

%Apply the scaling of the unit cell before differentiation
A = TPMS.tform.A(1:3,1:3);

xyz = inv(A)*([x y z]'-TPMS.tform.A(1:3,4));

% Create implicit functions for MC and GC using differential geometry
G = gradient(TPMS.u(xyz(1),xyz(2),xyz(3)),[x,y,z]);
G2 = [gradient(G(1),[x y z]) gradient(G(2),[x y z]) gradient(G(3),[x y z])];
K = det([G2' G; G' 0]);
H = -1/2*divergence(G/norm(G));
MCf = matlabFunction(H,'Vars',[x y z]);
GCf = matlabFunction(-K/(G'*G).^2,'Vars',[x y z]);

%Calculate the properties at the transformed points of interest
GC = GCf(points(:,1),points(:,2),points(:,3));
MC = MCf(points(:,1),points(:,2),points(:,3));
k1 = real(MC+sqrt(MC.^2-GC));
k2 = real(MC-sqrt(MC.^2-GC));
end

