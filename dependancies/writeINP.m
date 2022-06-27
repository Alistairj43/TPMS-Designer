function writeINP(elementStruct,nodeStruct,filename)
%% Function to write mesh data to INP
% Inputs:
% nodeStruct : A structure array containing the following fields: 
%       nodeStruct.N :An nx3 array of nodal coordinates
%       nodeStruct.N_ind :An nx1 array of the node indices (numbers)
%
% elementStruct : A a structure array of the form (for example):
%       elementStruct.E_type: '*ELEMENT, TYPE=C3D8R, ELSET=TPMSDesigner-Elements'  
%               String Specifying element type as the ABAQUS string
%       elementStruct.E: An mxl array of the nodal connectivity
%       elementStruct.E_ind: An mx1 array for the element indices
% filename : Name of the file to save to
% Outputs: saved to filename.inp file for intergration with ABAQUS (or similar)
%
% format of Inp temaplate file
% ------------------------------------------------------------------------
% Heading
% Parts
% Nodes
% Elements
% ------------------------------------------------------------------------
% Created by Alistair Jones, RMIT University 2022.



%Load from template INP file
fid=fopen('data/inp/template.inp');
T=textscan(fid,'%s','delimiter', '\n','Whitespace','');
T=T{1,1};
fclose(fid);

% Finding node and element fields in tempalte INP file
targets={'*NODE','*ELEMENT','**'}; %Text form targets
lineCount=1;
lineIndexTarget=zeros(size(targets));
i_target=1;
while 1
    l=T{lineCount};
    target=targets{i_target};
    if (strfind(l,target))
        lineIndexTarget(i_target)=lineCount;
        i_target=i_target+1;
    end    
    lineCount=lineCount+1;
    if nnz(lineIndexTarget)==(numel(targets))
        break
    end
end

%Creating new node text field
NODE_field=cell(size(nodeStruct.N,1),1);
for q=1:1:size(nodeStruct.N,1)
    nodeFieldLine=[sprintf('%8d,',nodeStruct.N_ind(q)) sprintf('% 10.6e, ',nodeStruct.N(q,:)) ];   
    NODE_field{q}=nodeFieldLine(1:end-2);    
end
 
%Creating new element text field
ELEMENT_field=cell(size(elementStruct.E,1),1);
for q=1:1:size(elementStruct.E,1)    
    elementFieldLine=sprintf('%8d,',[elementStruct.E_ind(q) elementStruct.E(q,:)]);   
    ELEMENT_field{q}=elementFieldLine(1:end-1);    
end

%% Generate the text
%Change element type line
T(lineIndexTarget(2))={elementStruct.E_type};

%Update node and element text fields
T=[T(1:lineIndexTarget(1)); NODE_field; T(lineIndexTarget(2)-2:lineIndexTarget(2)); ELEMENT_field; T(lineIndexTarget(3):end)];

%% Write to file
fid=fopen(filename,'wt');
for q=1:size(T,1)
    fprintf(fid,'%s\n',T{q});
end
fclose(fid);
end