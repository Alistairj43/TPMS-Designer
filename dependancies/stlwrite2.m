function stlwrite2(F,V,C)
cdata = [app.FV.currentProperty];
cLims = [min(cdata) max(cdata)];      % Transform height values
nCols = 256;  cMap = jet(nCols);      % onto an 16-bit colour map
fColsDbl = interp1(linspace(cLims(1),cLims(2),nCols),cMap,cdata);
if ~isempty(app.FVcap.faces)
    F = [app.FV.faces; app.FVcap.faces+length(app.FV.vertices)];
    V = [app.FV.vertices; app.FVcap.vertices];
    rgb = [uint16(fColsDbl*255); zeros(length(app.FVcap.faces),3)]; % Pass cols in 8bit (0-255) RGB triplets -Endcaps 0,0,0
else
    F = app.FV.faces;
    V = app.FV.vertices;
    rgb = uint16(fColsDbl*255);
end

switch event.Source
    case app.ColoredMeshplyMenu
        filename = "Mesh_Files/"+app.LabelEditField.Value+".ply";
        plywrite(filename,F,V,rgb);
    otherwise
        filename = "Mesh_Files/"+app.LabelEditField.Value+".stl";
        stlwrite(filename,F,V);
end
end