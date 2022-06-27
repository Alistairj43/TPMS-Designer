function h =plotVoxel(ax,data,property,opt)
%% Function plots a voxelated patch, coloured by "property"
% Inputs:   ax              - Handle to axis to plot on
%           data            - surfaceMesh object containing FV data
%           property        - string, the name of the property to be plotted
% Outputs:  h               - handle to the created patch graphical object
%
%   Credit: Alistair Jones, 2022, RMIT University

try
    pField = data.property.(extractBefore(property+' ',' '));
catch
    property = "solid";
    pField = data.property.solid;
end

voxels = data.property.solid(1:(data.res(1)-1),1:(data.res(2)-1),1:(data.res(3)-1));
cData = pField(1:(data.res(1)-1),1:(data.res(2)-1),1:(data.res(3)-1));

%detect the external voxels and faces
[FV, cData] = FindExternalVoxels(voxels,cData);
%plot only external faces of external voxels
if size(FV.vertices,1)==0
    h=[];
else
    FV.vertices = FV.vertices.*(data.upper-data.lower)./(data.res-1);
    h=patch(ax,'Faces',FV.faces,'Vertices',FV.vertices,'FaceVertexCData',cData,'FaceColor','flat');
    xlabel(ax,"X"); ylabel(ax,"Y"); zlabel(ax,"Z");
    c=colorbar(ax,'EastOutside'); c.Label.String = property; colormap(ax,"jet");
    rotate3d(ax,'on');
    
    if isfield(opt,'fancy')&&opt.fancy
        material(h,'shiny');
        view(ax,62,31);
        lighting(ax,'gouraud');
        ax.XTick = []; ax.YTick = []; ax.ZTick = [];
        ax.BoxStyle = 'full';
        camlight(ax);
        axis(ax,'manual','vis3d','equal','tight');
    end
end
end

function [FV, cData]=FindExternalVoxels(VoxelMat,cin)
% FindExternalVoxels scans VoxeLMat (a 3D matrix) and finds which voxels
% are external by checking if they have neighbors from all 6 sides
% (TOP,BOTTOM,FRONT,BACK,LEFT,RIGHT). After finding if the voxel is
% external it finds which faces are external by calling FindExternalFaces
% function and generating an FV structure for visualization
%initializing variables
FV.vertices=zeros(8*size(VoxelMat,1)*size(VoxelMat,2)*size(VoxelMat,3),3);
FV.faces=zeros(6*size(VoxelMat,1)*size(VoxelMat,2)*size(VoxelMat,3),4);
FaceIndex=1;
VertexIndex=0;
counter=1;
ExternalIndexes=zeros(size(VoxelMat,1)*size(VoxelMat,2)*size(VoxelMat,3),3);
voxel_size=[1 1 1];
for i=1:size(VoxelMat,1)
    for j=1:size(VoxelMat,2)
        for k=1:size(VoxelMat,3)
            if VoxelMat(i,j,k)==1
                if i==1 || j==1 || k==1 || i==size(VoxelMat,1) || j== size(VoxelMat,2) || k== size(VoxelMat,3)                 
                    ExternalIndexes(counter,1:3)=[i j k] ;
                    [FV,FaceIndex,VertexIndex]=FindExternalFaces(VoxelMat,ExternalIndexes,voxel_size,counter,FaceIndex,VertexIndex,FV);
                    cData(8*counter-7:8*counter,1) = cin(i,j,k)*ones(8,1);
                    counter=counter+1;
                else
                    if VoxelMat(i+1,j,k)==0 || VoxelMat(i-1,j,k)==0 || VoxelMat(i,j+1,k)==0 || VoxelMat(i,j-1,k)==0 || VoxelMat(i,j,k+1)==0 || VoxelMat(i,j,k-1)==0
                        ExternalIndexes(counter,1:3)=[i j k] ;
                        [FV,FaceIndex,VertexIndex]=FindExternalFaces(VoxelMat,ExternalIndexes,voxel_size,counter,FaceIndex,VertexIndex,FV);
                        cData(8*counter-7:8*counter,1) = cin(i,j,k)*ones(8,1);
                        counter=counter+1;
                    end
                end
            end
        end
    end 
end

FV.vertices=FV.vertices(any(FV.vertices,2),:);
FV.faces=FV.faces(any(FV.faces,2),:);
end

function [FV,FaceIndex,VertexIndex]=FindExternalFaces(VoxelMat,ExternalIndexes,voxel_size,i,FaceIndex,VertexIndex,FV)
faces=[1 2 3 4; 2 6 7 3 ; 6 5 8 7; 5 1 4 8; 4 3 7 8 ; 1 2 6 5];
FV.vertices(VertexIndex+1:VertexIndex+8,:)=[(ExternalIndexes(i,1)-1)+[0 voxel_size(1) voxel_size(1) 0 0 voxel_size(1) voxel_size(1) 0]; ...
    (ExternalIndexes(i,2)-1)+[0 0 0 0 voxel_size(2) voxel_size(2) voxel_size(2) voxel_size(2)]; ...
    (ExternalIndexes(i,3)-1)+[0 0 voxel_size(3) voxel_size(3) 0 0 voxel_size(3) voxel_size(3)]]';
if ExternalIndexes(i,2)~=1
    if VoxelMat(ExternalIndexes(i,1),ExternalIndexes(i,2)-1,ExternalIndexes(i,3))==0
        %No Front neighbor
        FV.faces(FaceIndex,:)=faces(1,:)+VertexIndex;
        FaceIndex=FaceIndex+1;
    end
else
    % Bounding Box Front
    FV.faces(FaceIndex,:)=faces(1,:)+VertexIndex;
    FaceIndex=FaceIndex+1;
end
if ExternalIndexes(i,1)~=size(VoxelMat,1)
    if VoxelMat(ExternalIndexes(i,1)+1,ExternalIndexes(i,2),ExternalIndexes(i,3))==0
        %No Right neighbor
        FV.faces(FaceIndex,:)=faces(2,:)+VertexIndex;
        FaceIndex=FaceIndex+1;
    end
else
    % Bounding Box Right
    FV.faces(FaceIndex,:)=faces(2,:)+VertexIndex;
    FaceIndex=FaceIndex+1;
end
if ExternalIndexes(i,2)~=size(VoxelMat,2)
    if VoxelMat(ExternalIndexes(i,1),ExternalIndexes(i,2)+1,ExternalIndexes(i,3))==0
        % No Back neighbor
        FV.faces(FaceIndex,:)=faces(3,:)+VertexIndex;
        FaceIndex=FaceIndex+1;
    end
else
    % Bounding Box Back
    FV.faces(FaceIndex,:)=faces(3,:)+VertexIndex;
    FaceIndex=FaceIndex+1;
end
if ExternalIndexes(i,1)~=1
    if VoxelMat(ExternalIndexes(i,1)-1,ExternalIndexes(i,2),ExternalIndexes(i,3))==0
        %No Left neighbor
        FV.faces(FaceIndex,:)=faces(4,:)+VertexIndex;
        FaceIndex=FaceIndex+1;
    end
else
    % Bounding Box Left
    FV.faces(FaceIndex,:)=faces(4,:)+VertexIndex;
    FaceIndex=FaceIndex+1;
end
if ExternalIndexes(i,3)~=size(VoxelMat,3)
    if VoxelMat(ExternalIndexes(i,1),ExternalIndexes(i,2),ExternalIndexes(i,3)+1)==0
        %No Top neighbor
        FV.faces(FaceIndex,:)=faces(5,:)+VertexIndex;
        FaceIndex=FaceIndex+1;
    end
else
    % Bounding Box Top
    FV.faces(FaceIndex,:)=faces(5,:)+VertexIndex;
    FaceIndex=FaceIndex+1;
end
if ExternalIndexes(i,3)~=1
    if VoxelMat(ExternalIndexes(i,1),ExternalIndexes(i,2),ExternalIndexes(i,3)-1)==0
        %No Bottom neighbor
        FV.faces(FaceIndex,:)=faces(6,:)+VertexIndex;
        FaceIndex=FaceIndex+1;
    end
else
    % Bounding Box Bottom
    FV.faces(FaceIndex,:)=faces(6,:)+VertexIndex;
    FaceIndex=FaceIndex+1;
end
VertexIndex=VertexIndex+8;
end