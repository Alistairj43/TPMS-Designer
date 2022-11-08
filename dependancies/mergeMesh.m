function [Fout,Vout,Cout]=mergeMesh(varargin)

% function [FT,VT,CT]=mergeMesh(Fc,Vc,Cc)
%------------------------------------------------------------------------

%% Parse inputs
switch nargin
    case 2
        Fc=varargin{1};
        Vc=varargin{2};        
        Cc=[];
    case 3
        Fc=varargin{1};
        Vc=varargin{2};
        Cc=varargin{3};
end

%%
numSets=numel(Fc); %Number of sets in the cell arrays

%Vertex data
numVertSets=cellfun(@(x) size(x,1),Vc); 
numDim=max(cellfun(@(x) size(x,2),Vc)); 
numVertTotal=sum(numVertSets); 
VT=zeros(numVertTotal,numDim); %Initialize

%Face data
numFaceSets=cellfun(@(x) size(x,1),Fc); 
numFaceVert=cellfun(@(x) size(x,2),Fc);
numFaceVertMax=max(numFaceVert);
numFaceVertMin=min(numFaceVert);
numFaceTotal=sum(numFaceSets); 
if numFaceVertMax==numFaceVertMin
   cellMode=0;
   FT=zeros(numFaceTotal,numFaceVertMax); %Initialize
else
    cellMode=1;
    FT=Fc; %Initialize as input cell
end

%Color data
if isempty(Cc) %If empty replace by set number
    Cc=Fc;
    for q=1:1:numSets
        Cc{q}=q*ones(size(Fc{q},1),1);
    end
end
numColorSets=cellfun(@(x) size(x,1),Cc);
numDimColor=max(cellfun(@(x) size(x,2),Cc));
numColorTotal=sum(numColorSets);

if cellMode        
    CT=Cc; %Initialize as input cell
else
    CT=zeros(numColorTotal,numDimColor); %Initialize
end

%% 

%Create index sets
indVert=[0; cumsum(numVertSets(:))];
indFace=[0; cumsum(numFaceSets(:))];
indColor=[0; cumsum(numColorSets(:))];

for q=1:1:numSets
    
    %Vertex data
    Vn=Vc{q}; %Current vertex set        
    if size(Vn,2)<numDim %Expand if required
        Vn(:,size(Vn,2)+1:numDim)=0; 
        warning('Not all vertex sets are of equal dimensionality (e.g. mixed 2D, 3D data). Zeros were added for added dimensions');
    end
    VT(indVert(q)+1:indVert(q+1),:)=Vn; %Appending Current vertex set
    
    %Face data
    Fn=Fc{q}+indVert(q); %Current face set plus index shift
    if cellMode
        FT{q}=Fn; %Appending Current face set
    else        
        FT(indFace(q)+1:indFace(q+1),:)=Fn; %Appending Current face set
    end    

    %Color Data
    if ~cellMode
        Cn=Cc{q}; %Current color set
        if size(Cn,2)<numDimColor %Expand if required
            Cn(:,size(Cn,2)+1:numDimColor)=0;
            warning('Not all color sets are of equal dimensionality (e.g. varying numbers of columns). Zeros were added for added dimensions');
        end
        CT(indColor(q)+1:indColor(q+1),:)=Cn; %Appending Current color set
    end
end

Fout=FT;
Vout=VT;
Cout=CT;

end