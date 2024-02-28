function writeFEmeshtoINP(FEmeshIN, fileName, elementType)
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

FEMesh.Elements = FEmeshIN.Elements';
FEMesh.Nodes = FEmeshIN.Nodes';
if(~strncmpi(fileName(end:-1:1),'pni.',4))
    fileName = [fileName,'.inp'];
end
while(exist(fileName,'file'))
    s = input('File name already exists, are you sure to overwrite it (1 for yes and 0 for no)?');
    switch s
        case {'yes', 1}
            break;
        otherwise
            fileName = input('Please input a new file name:','s');
    end
end
if(~strncmpi(fileName(end:-1:1),'pni.',4))
    fprintf('File name not with suffix ''.inp'', program will add it\n')
    fileName = [fileName,'.inp'];
end
fid = fopen(fileName,'wt');
fprintf(fid,'*HEADING\n');
fprintf(fid,'%s\n**\n',fileName);
fprintf(fid,'*Node\n');
for i = 1:1:size(FEMesh.Nodes,1)
    ndim = size(FEMesh.Nodes,2);
    fprintf(fid, ['%d',repmat(',%f', 1, ndim),'\n'], i, FEMesh.Nodes(i,:));
end
switch lower(elementType)
    case{'c2d3','c2d4','c2d6','c2d8'}
        eType = elementType;
        eType(2:3) = 'ps';
    otherwise
        eType = elementType;
end
fprintf(fid,'*Element, Type=%s\n', upper(eType));
for i = 1:1:size(FEMesh.Elements,1)
    nenode = size(FEMesh.Elements,2);
    fprintf(fid,['%d',repmat(',%d', 1, nenode),'\n'], i, FEMesh.Elements(i,:));
end
fclose(fid);
fprintf('Output Inp File Successfully!\n')
end