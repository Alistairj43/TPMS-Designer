function plotField(ax,data,zslice)
%% Function for adding a slice of the field
% Inputs:   N  (n x 2or3)    - Orientation vector
%     if (n x 2) use spherical cs values in radians [azimuth,inclination]
%     if (n x 3) use cartesian values[Nx,Ny,Nz]
%           Aw  (n x 1)       - Data Weighting
%           n  (int)         - number of subdivisons for icosphere mesh
% Outputs:  Vs (nv x 3)      - verticies for patch
%           Fs (nf x 3)      - faces for patch
%           Cs (nv x 1)      - per vertex colour data
% Designed for outputs to be used in matlab patch graphic,
%   ex: patch('Faces',Fs,'Vertices',Vs,'CData',Cs,'FaceColor','interp')
%
%   Credit: Alistair Jones, 12/11/2020, RMIT University
% Handle incorrect inputs
try
    pField = data.field.(extractBefore(data.p1,' '));
catch
    data.p1 = "U";
    pField = data.field.(data.p1);
end

z = range(data.field.zq);
[xq, yq, zq] = meshgrid(data.field.xq,data.field.yq,z*zslice/100);
xyslice = interp3(data.field.X,data.field.Y,data.field.Z,pField,xq,yq,zq);
surf(ax,xq,yq,zq,xyslice,'LineStyle','none');
[S,L] = bounds(pField,'all','omitnan');
c=colorbar; c.Label.String = data.p1; colormap(ax,gray);
caxis(ax,'manual'); caxis(ax,[S L]);
end