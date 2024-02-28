function writeFEmeshtoINP(FEmeshIN, filename, elementType)
%% Function to write mesh data to INP
% Inputs:
% femeshIN : Matlab PDE FEmesh, must have Nodes and Elements fields
%
% fileName : Name of the file to save to
% elementType : Type of elements
% Outputs: None - saves to filename.inp file for intergration with ABAQUS (or similar)
%
% format of .INP temaplate file
% ------------------------------------------------------------------------
% Heading
% Parts
% Nodes
% Elements
% ------------------------------------------------------------------------
% Original credit Wan Ji, 2021, 
% See: https://au.mathworks.com/matlabcentral/answers/1438724-how-to-export-femesh
% Adapted for use in TPMSDesigner by Alistair Jones, 2024.
arguments
    FEmeshIN;
    filename string = 'myFEMesh.inp';
    elementType string = 'c3d10';
end

[filepath,name,ext] = fileparts(filename);

FEMesh.Elements = FEmeshIN.Elements';
FEMesh.Nodes = FEmeshIN.Nodes';
fid = fopen(filename,'wt');
fprintf(fid,'*HEADING\n');
fprintf(fid,'%s\n**\n',[name '.inp']);
fprintf(fid,'*Node\n');
for i = 1:1:size(FEMesh.Nodes,1)
    ndim = size(FEMesh.Nodes,2);
    fprintf(fid, ['%d',repmat(',%f', 1, ndim),'\n'], i, FEMesh.Nodes(i,:));
end

fprintf(fid,'*Element, Type=%s\n', upper(elementType));
for i = 1:1:size(FEMesh.Elements,1)
    nenode = size(FEMesh.Elements,2);
    fprintf(fid,['%d',repmat(',%d', 1, nenode),'\n'], i, FEMesh.Elements(i,:));
end
fclose(fid);
fprintf('Output Inp File Successfully!\n')
end