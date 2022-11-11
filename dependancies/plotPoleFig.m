function h = plotPoleFig(FV,n,opts,ax)
%% Function for plotting a heatmap of orientation, mapped onto a unit sphere surface
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
arguments
    FV;
    n = 3;
    opts = [];
    ax = [];
end

if isempty(ax)
    figure; ax = gca; view(3);
end

if isempty(FV.faces)
    h = [];
else
    N = [FV.Fproperty.Nx FV.Fproperty.Ny FV.Fproperty.Nz];
    Aw = FV.Fproperty.Farea;

    totalA = sum(Aw,'all','omitnan'); %Total per face weighting
    [FV] = icosphere(n); %Generate an icosphere with n subdivisions
    ids = dsearchn(FV.vertices,N); %Find which node of the icosphere is closest to the unit normal vector
    Cs = 4*pi*length(FV.vertices)*accumarray(ids,Aw,[length(FV.vertices) 1],@sum,0)/totalA+0.001; %Complete cumulative sum normalise by area of sphere and surface

    %surf(ax,Fc(:,1),Fc(:,2),Fc(:,3),Cs,'FaceColor','interp','EdgeColor','none')
    h = patch(ax,'Faces',FV.faces,'Vertices',FV.vertices,'CData',Cs,'FaceColor','interp','EdgeColor','none');
    axis(ax,'equal','vis3d'); view(ax,3);

    c=colorbar(ax); c.Label.String = 'Relative Intensity'; colormap(ax,"jet");
    %caxis(ax,prctile(Cs,[0.2 99.8]));
    axis(ax,'manual','vis3d','equal','tight');
    xlabel(ax,"X"); ylabel(ax,"Y"); zlabel(ax,"Z");
    ax.XTick = []; ax.YTick = []; ax.ZTick = [];
    ax.BoxStyle = 'full';
    rotate3d(ax,'on');
    view(ax,62,31);
end
end