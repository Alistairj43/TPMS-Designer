function plotMesh(ax,data,property)
%% Function returns points for a heatmap mapped onto a unit sphere surface
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
opts.LineStyle='none';

try
    pField = data.property.(extractBefore(property+' ',' '));
catch
    p1 = "inclination angle (^o deg)";
    pField = data.property.inclination;
end

if ~isempty(data.FVcap.faces)
    hold(ax,"on"); patch(ax,'Faces',data.facesC,'Vertices',data.verticesC,'FaceColor',[0.25,0.25,0.25],...
        'FaceAlpha',0.9,opts);
end

if ~isempty(data.FV.faces)
    opts.FaceColor = 'flat';
    patch(ax,'Faces',data.faces,'Vertices',data.vertices,'FaceVertexCData',pField,opts);
end

xlabel(ax,"X (mm)"); ylabel(ax,"Y (mm)"); zlabel(ax,"Z - Build Direction (mm)");
c=colorbar; c.Label.String = p1; colormap(ax,"jet");
end