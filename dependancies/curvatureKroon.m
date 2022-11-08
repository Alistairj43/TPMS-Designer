function [GC, MC, k1, k2, ABS, RMS]=curvatureKroon(F,V,usethird)
% This function calculates the principal curvature directions and values
% of a triangulated mesh. 
% [Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2,angle]=patchcurvature(F,V,usethird)
%
% inputs,
%   FV : A triangulated mesh (see Patch)
%   usethird : Use third order neighbour vertices for the curvature
%              fit, making it smoother but less local. true/ false (default)
% Outputs:
%   GC - gaussian curvature (from determinant of the shape matrix)
%   MC - mean curvature (from divergance of unit normal)
%   k1 - maximum principle curvature (from GC and MC)
%   k2 - minimum principle curvature (from GC and MC)
%   ABS - absoloute curvature (from k1 and k2)
%   RMS - rms curvature (from GC and MC)
% Designed for analysis of patch objects
%
% Modified by Alistair Jones, 30/09/2021, RMIT University,
% Original Function is written by D.Kroon University of Twente (August 2011)  

% Check inputs
if(nargin<3), usethird=true; end
% Number of vertices
nv=size(V,1);
% Calculate vertices normals
N=patchnormals(F,V);

% Calculate Rotation matrices for the normals
M= zeros(3,3,nv);
Minv= zeros(3,3,nv);
for i=1:nv 
    [M(:,:,i),Minv(:,:,i)]=VectorRotationMatrix(N(i,:));
end

% Get neighbours of all vertices
Ne=vertex_neighbours(F,V);

% Loop through all vertices
k1=zeros(nv,1);
k2=zeros(nv,1);
Dir1=zeros(nv,3);
Dir2=zeros(nv,3);

for i=1:nv
   % Get first and second ring neighbours.
   if(~usethird)
       Nce=unique([Ne{Ne{i}}]);
   else
       % Get first, second and third ring neighbours
       Nce=unique([Ne{[Ne{Ne{i}}]}]);
   end
   
   Ve=V(Nce,:);

   % Rotate to make normal [-1 0 0]
   We=Ve*Minv(:,:,i);
   f=We(:,1); x=We(:,2); y=We(:,3); 
   
   % Fit patch
   % f(x,y) = ax^2 + by^2 + cxy + dx + ey + f
   FM=[x(:).^2 y(:).^2 x(:).*y(:) x(:) y(:) ones(numel(x),1)];
   abcdef=FM\f(:);
   a=abcdef(1); b=abcdef(2); c=abcdef(3);
   
   % Make Hessian matrix 
   % H =  [2*a c;c 2*b];
   Dxx = 2*a; Dxy=c; Dyy=2*b;
   
   [k1(i),k2(i),I1,I2]=eig2(Dxx,Dxy,Dyy);

   dir1=[0 I1(1) I1(2)]*M(:,:,i); 
   dir2=[0 I2(1) I2(2)]*M(:,:,i);
   Dir1(i,:)=dir1/sqrt(dir1(1)^2+dir1(2)^2+dir1(3)^2);
   Dir2(i,:)=dir2/sqrt(dir2(1)^2+dir2(2)^2+dir2(3)^2);
end

MC=(k1+k2)/2;
GC=k1.*k2;
ABS = abs(k1)+abs(k2);
RMS = real(sqrt(4*MC.^2-2*GC));
end 


function [Lambda1,Lambda2,I1,I2]=eig2(Dxx,Dxy,Dyy)
% Compute the eigenvectors 
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
v2x = 2*Dxy; v2y = Dyy - Dxx + tmp;

% Normalize
mag = sqrt(v2x.^2 + v2y.^2); i = (mag ~= 0);
v2x(i) = v2x(i)./mag(i);
v2y(i) = v2y(i)./mag(i);

% The eigenvectors are orthogonal
v1x = -v2y; v1y = v2x;

% Compute the eigenvalues
mu1 = (0.5*(Dxx + Dyy + tmp));
mu2 = (0.5*(Dxx + Dyy - tmp));

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
if(abs(mu1)<abs(mu2))
    Lambda1=mu1;
    Lambda2=mu2;
    I2=[v1x v1y];
    I1=[v2x v2y];
else
    Lambda1=mu2;
    Lambda2=mu1;
    I2=[v2x v2y];
    I1=[v1x v1y];
end
end

function N=patchnormals(F,V)
% This function PATCHNORMALS calculates the normals of a triangulated
% mesh. PATCHNORMALS calls the patchnormal_double.c mex function which 
% first calculates the normals of all faces, and after that calculates 
% the vertice normals from the face normals weighted by the angles 
% of the faces.
[Nx,Ny,Nz]=patchnormals_double(double(F(:,1)),double(F(:,2)),double(F(:,3)),double(V(:,1)),double(V(:,2)),double(V(:,3)));
N=zeros(length(Nx),3);
N(:,1)=Nx; N(:,2)=Ny; N(:,3)=Nz;
end

function [Nx,Ny,Nz]=patchnormals_double(Fa,Fb,Fc,Vx,Vy,Vz)
%
%  [Nx,Ny,Nz]=patchnormals_double(Fa,Fb,Fc,Vx,Vy,Vz)
%
FV.vertices=zeros(length(Vx),3);
FV.vertices(:,1)=Vx;
FV.vertices(:,2)=Vy;
FV.vertices(:,3)=Vz;

% Get all edge vectors
e1=FV.vertices(Fa,:)-FV.vertices(Fb,:);
e2=FV.vertices(Fb,:)-FV.vertices(Fc,:);
e3=FV.vertices(Fc,:)-FV.vertices(Fa,:);

% Normalize edge vectors
e1_norm=e1./repmat(sqrt(e1(:,1).^2+e1(:,2).^2+e1(:,3).^2),1,3); 
e2_norm=e2./repmat(sqrt(e2(:,1).^2+e2(:,2).^2+e2(:,3).^2),1,3); 
e3_norm=e3./repmat(sqrt(e3(:,1).^2+e3(:,2).^2+e3(:,3).^2),1,3);

% Calculate Angle of face seen from vertices
Angle =  [acos(dot(e1_norm',-e3_norm'));acos(dot(e2_norm',-e1_norm'));acos(dot(e3_norm',-e2_norm'))]';

% Calculate normal of face
 Normal=cross(e1,e3);

% Calculate Vertice Normals 
VerticeNormals=zeros([size(FV.vertices,1) 3]);
for i=1:size(Fa,1)
    VerticeNormals(Fa(i),:)=VerticeNormals(Fa(i),:)+Normal(i,:)*Angle(i,1);
    VerticeNormals(Fb(i),:)=VerticeNormals(Fb(i),:)+Normal(i,:)*Angle(i,2);
    VerticeNormals(Fc(i),:)=VerticeNormals(Fc(i),:)+Normal(i,:)*Angle(i,3);
end

V_norm=sqrt(VerticeNormals(:,1).^2+VerticeNormals(:,2).^2+VerticeNormals(:,3).^2)+eps;
VerticeNormals=VerticeNormals./repmat(V_norm,1,3);
Nx=VerticeNormals(:,1);
Ny=VerticeNormals(:,2);
Nz=VerticeNormals(:,3);
end

function [M,Minv]=VectorRotationMatrix(v)
% [M,Minv]=VectorRotationMatrix(v,k)
v=(v(:)')/sqrt(sum(v.^2));
k=rand(1,3);
l = [k(2).*v(3)-k(3).*v(2), k(3).*v(1)-k(1).*v(3), k(1).*v(2)-k(2).*v(1)]; l=l/sqrt(sum(l.^2));
k = [l(2).*v(3)-l(3).*v(2), l(3).*v(1)-l(1).*v(3), l(1).*v(2)-l(2).*v(1)]; k=k/sqrt(sum(k.^2));
Minv=[v(:) l(:) k(:)];
M=inv(Minv);
end

function Ne=vertex_neighbours(F,V)
% This function VERTEX_NEIGHBOURS will search in a face list for all 
% the neigbours of each vertex.
%
% Ne=vertex_neighbours(FV)
%
Ne=vertex_neighbours_double(F(:,1),F(:,2),F(:,3),V(:,1),V(:,2),V(:,3));
end

function Ne=vertex_neighbours_double(Fa,Fb,Fc,Vx,Vy,Vz)

F=[Fa Fb Fc];
V=[Vx Vy Vz];

% Neighbourh cell array 
Ne=cell(1,size(V,1));

% Loop through all faces
for i=1:length(F)
    % Add the neighbors of each vertice of a face
    % to his neighbors list.
    Ne{F(i,1)}=[Ne{F(i,1)} [F(i,2) F(i,3)]];
    Ne{F(i,2)}=[Ne{F(i,2)} [F(i,3) F(i,1)]];
    Ne{F(i,3)}=[Ne{F(i,3)} [F(i,1) F(i,2)]];
end

% Loop through all neighbor arrays and sort them (Rotation same as faces)
for i=1:size(V,1)
    Pneighf=Ne{i};
    if(isempty(Pneighf))
        Pneig=[];
    else
        start=1;
        for index1=1:2:length(Pneighf)
            found=false;
            for index2=2:2:length(Pneighf)
                if(Pneighf(index1)==Pneighf(index2))
                    found=true; break
                end
            end
            if(~found)
                start=index1; break
            end
        end
        Pneig=[];
        Pneig(1)=Pneighf(start);
        Pneig(2)=Pneighf(start+1);
        
        % Add the neighbours with respect to original rotation
        for j=2+double(found):(length(Pneighf)/2)
            found = false;
            for index=1:2:length(Pneighf)
                if(Pneighf(index)==Pneig(end))
                    if(sum(Pneig==Pneighf(index+1))==0)
                        found =true;
                        Pneig=[Pneig Pneighf(index+1)];
                    end
                end
            end
            if(~found) % This only happens with weird edge vertices
                for index=1:2:length(Pneighf)
                    if(sum(Pneig==Pneighf(index))==0)
                        Pneig=[Pneig Pneighf(index)];
                        if(sum(Pneig==Pneighf(index+1))==0)
                            Pneig=[Pneig Pneighf(index+1)];
                        end
                    end
                end
            end
        end
        % Add forgotten neigbours
        if(length(Pneig)<length(Pneighf))
            for j=1:length(Pneighf)
                if(sum(Pneig==Pneighf(j))==0)
                    Pneig=[Pneig Pneighf(j)];
                end
            end
        end
    end
    Ne{i}=Pneig;
end
end
