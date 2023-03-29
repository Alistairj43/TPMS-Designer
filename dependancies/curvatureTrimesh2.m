function [kout1,kout2,kDir1,kDir2,wfp]=curvatureTrimesh2(FV,VertexNormals,FaceNormals,Avertex,Acorner,up,vp)
%% Summary
%recives a list of vertices and faces in FV structure
%and the normal at each vertex and calculates the second fundemental
%matrix and the curvature using least squares
%INPUT:
%FV - face-vertex data structure containing a list of vertices and a list of faces
%VertexNormals - nX3 matrix (n=number of vertices) containing the normal at each vertex
%FaceNormals - mX3 matrix (m = number of faces) containing the normal of each face
%OUTPUT:
%FaceSFM - an mX1 cell matrix (m = number of faces) second fundemental
%VertexSFM - an nX1 cell matrix (n = number of vertices) second fundemental
%wfp - corner voronoi weights
%% Code
%matrix of each face at each cell
FaceSFM=cell(size(FV.faces,1),1);
VertexSFM=cell(size(FV.vertices,1),1);
% [FaceSFM{1:end,1}]=deal(zeros(2,2));
% [VertexSFM{1:end,1}]=deal(zeros(2,2));
FaceSFM(1:end,1)={zeros(2,2)};
VertexSFM(1:end,1)={zeros(2,2)}; 
% Get all edge vectors
e0=FV.vertices(FV.faces(:,3),:)-FV.vertices(FV.faces(:,2),:);
e1=FV.vertices(FV.faces(:,1),:)-FV.vertices(FV.faces(:,3),:);
e2=FV.vertices(FV.faces(:,2),:)-FV.vertices(FV.faces(:,1),:);
% Normalize edge vectors
e0_norm=normalize(e0,2,'norm');

wfp=zeros(size(FV.faces,1),3);

for i=1:size(FV.faces,1)
    %Calculate Curvature Per Face
    %set face coordinate frame
    nf=FaceNormals(i,:);
    t=e0_norm(i,:)';
    B=cross(nf,t)';
    B= B/norm(B);
    %extract relevant normals in face vertices
    n0=VertexNormals(FV.faces(i,1),:);
    n1=VertexNormals(FV.faces(i,2),:);
    n2=VertexNormals(FV.faces(i,3),:);
    %solve least squares problem of th form Ax=b
    A=[e0(i,:)*t e0(i,:)*B 0;
        0 e0(i,:)*t e0(i,:)*B;
        e1(i,:)*t e1(i,:)*B 0;
        0 e1(i,:)*t e1(i,:)*B;
        e2(i,:)*t e2(i,:)*B 0;
        0 e2(i,:)*t e2(i,:)*B];
    b=[(n2-n1)*t;(n2-n1)*B;(n0-n2)*t;(n0-n2)*B;(n1-n0)*t;(n1-n0)*B];
    %[LA,DA] = ldl(A'*A);
    % bA=A'*b;
    %  x = LA'\(DA\(LA\bA));
    x=A\b;
    
    FaceSFM{i,1}=[x(1),x(2);x(2) x(3)];
    %Kn(i)=[1 0]*FaceSFM{i,1}*[1;0];
    
    %Calculate Curvature Per Vertex
    %calculate voronoi weights
    wfp(i,1)=Acorner(i,1)/Avertex(FV.faces(i,1));
    wfp(i,2)=Acorner(i,2)/Avertex(FV.faces(i,2));
    wfp(i,3)=Acorner(i,3)/Avertex(FV.faces(i,3));
    %Calculate new coordinate system and project the tensor
    for j=1:3
        [new_ku,new_kuv,new_kv]=ProjectCurvatureTensor(t,B,nf,x(1),x(2),x(3),up(FV.faces(i,j),:),vp(FV.faces(i,j),:));
        VertexSFM{FV.faces(i,j),1}= VertexSFM{FV.faces(i,j),1}+wfp(i,j)*[new_ku new_kuv;new_kuv new_kv];
    end
end


kout1=zeros(size(FV.vertices,1),1);
kout2=zeros(size(FV.vertices,1),1);
kDir1=zeros(size(FV.vertices,1),3);
kDir2=zeros(size(FV.vertices,1),3);
%{
%solve using matlab eigen value and eigenvector solver
for i=1:size(FV.vertices)
    [V,D]=eig(VertexSFM{i,1});
    % PrincipalCurvatures(:,i)=eig(VertexSFM{i,1});
    [maxk indexk]=max(abs([D(1,1),D(2,2)]));
    if indexk==1
        PrincipalCurvatures(1,i)=D(1,1);
        PrincipalCurvatures(2,i)=D(2,2);
        PrincipalDir1(i,:)=(V(1,1)*up(i,:)+V(2,1)*vp(i,:));
        PrincipalDir2(i,:)=(V(1,2)*up(i,:)+V(2,2)*vp(i,:));
    else
        PrincipalCurvatures(2,i)=D(1,1);
        PrincipalCurvatures(1,i)=D(2,2);
        PrincipalDir2(i,:)=(V(1,1)*up(i,:)+V(2,1)*vp(i,:));
        PrincipalDir1(i,:)=(V(1,2)*up(i,:)+V(2,2)*vp(i,:));
    end
end
%}
%
for i=1:size(FV.vertices)
%This is taken from trimsh2 - Szymon Rusinkiewicz implementation. 
% It also considers direction
np=cross(up(i,:),vp(i,:));
[r_old_u, r_old_v]=RotateCoordinateSystem(up(i,:), vp(i,:), np);
ku=VertexSFM{i,1}(1,1);
kuv=VertexSFM{i,1}(1,2);
kv=VertexSFM{i,1}(2,2);
c = 1;
s = 0;
tt = 0;
if kuv ~= 0
    %Jacobi rotation to diagonalize
    h = 0.5 * (kv - ku) / kuv;
    if  h < 0
        tt =1 / (h - sqrt(1 + h*h)) ;
    else
        tt =1/ (h + sqrt(1 + h*h));
    end
    c = 1 / sqrt(1+ tt*tt);
    s = tt * c;
end

k1 = ku - tt * kuv;
k2 = kv + tt * kuv;

if (abs(k1) >= abs(k2))
    kDir1(i,:) = c*r_old_u - s*r_old_v;
else
    temp=k1;
    k1=k2;
    k2=temp;
    kDir1(i,:) = s*r_old_u + c*r_old_v;
end
kDir2(i,:) = cross(np , kDir1(i,:));
kout1(i)=k1;
kout2(i)=k2;
end

end

function [new_ku,new_kuv,new_kv]=ProjectCurvatureTensor(uf,vf,nf,old_ku,old_kuv,old_kv,up,vp)
%{ Summary: ProjectCurvatureTensor performs a projection
%of the tensor variables to the vertexcoordinate system
%INPUT:
%uf,vf - face coordinate system
%old_ku,old_kuv,old_kv - face curvature tensor variables
%up,vp - vertex cordinate system
%OUTPUT:
%new_ku,new_kuv,new_kv - vertex curvature tensor variabels
%}
[r_new_u,r_new_v]=RotateCoordinateSystem(up,vp,nf);
OldTensor=[old_ku old_kuv; old_kuv old_kv];
u1=r_new_u*uf;
v1=r_new_u*vf;
u2=r_new_v*uf;
v2=r_new_v*vf;
new_ku=[u1 v1]*OldTensor*[u1;v1];
new_kuv=[u1 v1]*OldTensor*[u2;v2];
new_kv=[u2 v2]*OldTensor*[u2;v2];
%{
new_ku=old_ku*u1*u1+2*old_kuv*(u1*v1)+old_kv*v1*v1;
new_kuv=old_ku*u1*u2+old_kuv*(u1*v2+u2*v1)+old_kv*v1*v2;
new_kv=old_ku*u2*u2+2*old_kuv*(u2*v2)+old_kv*v2*v2;
%}
end

function [r_new_u,r_new_v]=RotateCoordinateSystem(up,vp,nf)
%{Summary: RotateCoordinateSystem performs the rotation of the vectors up and vp
%to the plane defined by nf
%INPUT:
%up,vp- vectors to be rotated (vertex coordinate system)
%nf - face normal
%OUTPUT:
%r_new_u,r_new_v - rotated coordinate system
%}
r_new_u=up;
r_new_v=vp;
np=cross(up,vp);
np=np/norm(np);
ndot=nf*np';
if ndot<=-1
    r_new_u=-r_new_u;
    r_new_v=-r_new_v;  
    return;
end
perp=nf-ndot*np;
dperp=(np+nf)/(1+ndot);
r_new_u=r_new_u-dperp*(perp*r_new_u');
r_new_v=r_new_v-dperp*(perp*r_new_v');
end
