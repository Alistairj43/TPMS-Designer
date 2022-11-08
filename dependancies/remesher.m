function [vnew, fnew, meanedge, stdev]=remesher(V, F, edgelength, iterations)
% INPUT
% -V: vertices of mesh; n*3 array of xyz coordinates
% -F: faces of  mesh; n*3 array
% -edgelength:  edgelength aimed for
% -iterations: nr of iterations, 5 to 10 suffice
% 
% OUTPUT
% -vnew: remeshed vertices list: n*3 array of xyz coordinates
% -fnew: remeshed list of faces of  mesh; n*3 array
% -meanedge: average edgelenth obtained
% -stdev: stdev of edgelengths
% EXAMPLE
%     load testdatafemur
%     [vnew, fnew, meanedge, stdev]=remesher(vnew, fnew, 3, 5);
%clean up patch
[vnew, fnew]=cleanpatch(V, F);
voriginal=vnew;
foriginal=fnew;
[vnew,fnew] = subdividelarge( vnew, fnew,edgelength,voriginal,foriginal );
voriginal=vnew;
foriginal=fnew;
for i=1:iterations
    
[vnew, fnew,temp] = edgecollaps( vnew, fnew, edgelength ,voriginal,foriginal);
[vnew, fnew, temp] = removebadtriangles( vnew, fnew,voriginal,foriginal);
[vnew,fnew] = subdividelarge( vnew, fnew,0 ,voriginal,foriginal);
[vnew, fnew, temp] = removebadtriangles( vnew, fnew,voriginal,foriginal);
disp(['Iteration:' num2str(i) '  Output mesh: ' num2str(size(fnew,1)) ' triangles, ' ... 
    num2str(size(vnew,1))  ' vertices.']);
end
meanedge=temp(:,1);
stdev=temp(:,2);
end

function [vnew, fnew,temp] = edgecollaps( vnew, fnew, sizzz,voriginal,foriginal )
[vnew, fnew]=cleanpatch(vnew, fnew);
cases=vnew;
sizzz=sizzz/sqrt(3);
while size(cases,1)>0
fk1 = fnew(:,1);
fk2 = fnew(:,2);
fk3 = fnew(:,3);
numfaces = (1:size(fnew,1))';
e1=sqrt(sum((vnew(fk1,:)-vnew(fk2,:)).^2,2));
e2=sqrt(sum((vnew(fk1,:)-vnew(fk3,:)).^2,2));
e3=sqrt(sum((vnew(fk2,:)-vnew(fk3,:)).^2,2));
temp(:,1)=mean(e1);
temp(:,2)=std(e1);
ed1=sort([fk1 fk2 ]')';
ed2=sort([fk1 fk3 ]')';
ed3=sort([fk2 fk3 ]')';
e1=[e1 numfaces ed1 ];
e2=[e2 numfaces ed2 ];
e3=[e3 numfaces ed3 ];
e=[e1 ; e2 ; e3];
e=e(e(:,1)<sizzz,:);
e=sortrows(e,1);
[etemp,ia,ic]=unique(e(:,3),'rows','stable');
e=e(ia,:);
[etemp,ia,ic]=unique(e(:,4),'rows','stable');
e=e(ia,:);
[test,ia,ic]=unique(e(:,2),'rows','stable');
e=e(ia,:);
ind1=(1:2:(2*size(ia,1)-1))';
ind2=(2:2:(2*size(ia,1)))';
ind3=(1:2*size(ia,1))';
test1=ones(2*size(ia,1),1);
test1(ind1)=e(:,3);
test1(ind2)=e(:,4);
test1(:,2)=ones;
test1(ind1,2)=(1:size(ia))';
test1(ind2,2)=(1:size(ia))';
[etemp1,ia,ic]=unique(test1(:,1),'stable');
test1=(test1(ia,:));
test1(:,3)=ones;
test1(2:end,3)=test1(1:end-1,2);
test1(:,4)=test1(:,3)-test1(:,2);
indicesseries= test1(test1(:,4)==0,2);
indicesseries=unique(indicesseries,'stable');
e=e(indicesseries,:);
% doubles= find(ismember(e(:,4),e(:,3)));
% e=removerows(e,doubles);
cases=e(:,3:4);
 
averages=(vnew(cases(:,1),:)+vnew(cases(:,2),:)).*0.5;
vnew(cases(:,1),:)=averages;
vnew(cases(:,2),:)=averages;
    
[vnew, fnew]=cleanpatch(vnew, fnew);
[vnew]=project(vnew, fnew,voriginal,foriginal);
end
end

function [projections]=project(vS,fS,vT,fT)
TRS = triangulation(fS,vS); 
normalsS=vertexNormal(TRS);
[IDXsource,Dsource]=knnsearch(vT,vS);
vector_s_to_t=vT(IDXsource,:)-vS;
projections=vS+[(sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,1) (sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,2) (sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,3)];
end

function [vnew, fnew, temp] = removebadtriangles( vnew, fnew ,voriginal,foriginal)
cases=vnew;
fk1 = fnew(:,1);
fk2 = fnew(:,2);
fk3 = fnew(:,3);
numfaces = (1:size(fnew,1))';
e1=sqrt(sum((vnew(fk1,:)-vnew(fk2,:)).^2,2));
e2=sqrt(sum((vnew(fk1,:)-vnew(fk3,:)).^2,2));
e3=sqrt(sum((vnew(fk2,:)-vnew(fk3,:)).^2,2));
%define area by Heron formula
s=(e1+e2+e3).*0.5;
area=sqrt(s.*(s-e1).*(s-e2).*(s-e3));
quality=[e1 e2 e3];
M = max(quality,[],2);
qualitycheck=area./M;
sizzz=mean(qualitycheck)-1.65*std(qualitycheck);
while size(cases,1)>0
fk1 = fnew(:,1);
fk2 = fnew(:,2);
fk3 = fnew(:,3);
numfaces = (1:size(fnew,1))';
e1=sqrt(sum((vnew(fk1,:)-vnew(fk2,:)).^2,2));
e2=sqrt(sum((vnew(fk1,:)-vnew(fk3,:)).^2,2));
e3=sqrt(sum((vnew(fk2,:)-vnew(fk3,:)).^2,2));
temp(:,1)=mean(e1);
temp(:,2)=std(e1);
%define area by Heron formula
s=(e1+e2+e3).*0.5;
area=sqrt(s.*(s-e1).*(s-e2).*(s-e3));
quality=[e1 e2 e3];
M = max(quality,[],2);
qualitycheck=area./M;
[m,idx]=min(quality,[],2);
score1=find(ismember(idx,[1]));
score2=find(ismember(idx,[2]));
score3=find(ismember(idx,[3]));
ed1=sort([fk1 fk2 ]')';
ed2=sort([fk1 fk3 ]')';
ed3=sort([fk2 fk3 ]')';
e1=[qualitycheck numfaces ed1 ones(size(e1,1),1)];
e2=[qualitycheck numfaces ed2 ones(size(e2,1),1)];
e3=[qualitycheck numfaces ed3 ones(size(e2,1),1)];
 e1(score1,5)=0;
 e2(score2,5)=0;
 e3(score3,5)=0;
 
e=[e1 ; e2 ; e3];
e=e(e(:,5)==0,1:4);
e=e(e(:,1)<sizzz,:);
e=sortrows(e,1);
[etemp,ia,ic]=unique(e(:,3),'rows','stable');
e=e(ia,:);
[etemp,ia,ic]=unique(e(:,4),'rows','stable');
e=e(ia,:);
[test,ia,ic]=unique(e(:,2),'rows','stable');
e=e(ia,:);
ind1=(1:2:(2*size(ia,1)-1))';
ind2=(2:2:(2*size(ia,1)))';
ind3=(1:2*size(ia,1))';
test1=ones(2*size(ia,1),1);
test1(ind1)=e(:,3);
test1(ind2)=e(:,4);
test1(:,2)=ones;
test1(ind1,2)=(1:size(ia))';
test1(ind2,2)=(1:size(ia))';
[etemp1,ia,ic]=unique(test1(:,1),'stable');
test1=(test1(ia,:));
test1(:,3)=ones;
test1(2:end,3)=test1(1:end-1,2);
test1(:,4)=test1(:,3)-test1(:,2);
indicesseries= test1(test1(:,4)==0,2);
indicesseries=unique(indicesseries,'stable');
e=e(indicesseries,:);
cases=e(:,3:4);
 
averages=(vnew(cases(:,1),:)+vnew(cases(:,2),:)).*0.5;
vnew(cases(:,1),:)=averages;
vnew(cases(:,2),:)=averages;
    
[vnew, fnew]=cleanpatch(vnew, fnew);
[vnew]=project(vnew, fnew,voriginal,foriginal);
end
end

function [vnew,fnew] = subdividelarge( vnew, fnew,flag,voriginal,foriginal )
vn=[vnew];
if flag==0
fk1 = fnew(:,1);
fk2 = fnew(:,2);
fk3 = fnew(:,3);
numfaces = (1:size(fnew,1))';
e1=sqrt(sum((vnew(fk1,:)-vnew(fk2,:)).^2,2));
e2=sqrt(sum((vnew(fk1,:)-vnew(fk3,:)).^2,2));
e3=sqrt(sum((vnew(fk2,:)-vnew(fk3,:)).^2,2));
sizzz=mean([e1;e2;e3])+1.96*std([e1;e2;e3]);
else
    
    sizzz=flag;
end
   
while size(vn,1)>1
fk1 = fnew(:,1);
fk2 = fnew(:,2);
fk3 = fnew(:,3);
numfaces = (1:size(fnew,1))';
e1=sqrt(sum((vnew(fk1,:)-vnew(fk2,:)).^2,2));
e2=sqrt(sum((vnew(fk1,:)-vnew(fk3,:)).^2,2));
e3=sqrt(sum((vnew(fk2,:)-vnew(fk3,:)).^2,2));
ed1=sort([fk1 fk2 ]')';
ed2=sort([fk1 fk3 ]')';
ed3=sort([fk2 fk3 ]')';
e1=[e1 numfaces ed1 fk3];
e2=[e2 numfaces ed2 fk2];
e3=[e3 numfaces ed3 fk1];
e=[e1 ; e2 ; e3];
e=e(e(:,1)>sizzz,:);
%single edges
e=sortrows(e,-1);
[etemp1,ia,ic]=unique(e(:,3:4),'rows','stable');
esingle=e(ia,:);
%dubbles
edouble=removerows(e,ia);
[C,ia,ib] = intersect(esingle(:,3:4),edouble(:,3:4),'rows','stable');
% newseries=[esingle(ia,:) edouble(ib,:)];
newseries=[esingle(ia,:) edouble(ib,:)];
newseries=sortrows(newseries,-1);
ind1=(1:2:(2*size(ia,1)-1))';
ind2=(2:2:(2*size(ia,1)))';
ind3=(1:2*size(ia,1))';
test1=ones(2*size(ia,1),1);
test1(ind1)=newseries(:,2);
test1(ind2)=newseries(:,7);
test1(:,2)=ones;
test1(ind1,2)=(1:size(ia))';
test1(ind2,2)=(1:size(ia))';
[etemp1,ia,ic]=unique(test1(:,1),'stable');
test1=(test1(ia,:));
test1(:,3)=ones;
test1(2:end,3)=test1(1:end-1,2);
test1(:,4)=test1(:,3)-test1(:,2);
indicesseries= test1(test1(:,4)==0,2);
indicesseries=unique(indicesseries,'stable');
newseries=newseries(indicesseries,:);
vn=(vnew(newseries(:,3),:)+vnew(newseries(:,4),:)).*0.5;
sizevn=size(vn,1);
indices=size(vnew,1)+(1:sizevn)';
e=[];
e=[horzcat(newseries(:,1:5),indices);horzcat(newseries(:,6:10),indices)];
faces=fnew(e(:,2),:);
n1=vnew(faces(:,1),:)-vnew(faces(:,2),:);
n3=vnew(faces(:,1),:)-vnew(faces(:,3),:);
Normals=cross(n1,n3);
Distance=sqrt(sum((Normals.^2),2));
Normalsoriginal=horzcat(Normals(:,1)./Distance,Normals(:,2)./Distance,Normals(:,3)./Distance);
Normalsoriginal=[Normalsoriginal;Normalsoriginal];
%define new vertices
vnew=[vnew ;vn];
sizevn=size(vn,1);
%define new faces
f1=[e(:,3) e(:,5) e(:,6)];
f2=[e(:,4) e(:,5) e(:,6)];
f=[f1 ; f2];
%correct normals
nn1=vnew(f(:,1),:)-vnew(f(:,2),:);
nn3=vnew(f(:,1),:)-vnew(f(:,3),:);
Normals=cross(nn1,nn3);
Distance=sqrt(sum((Normals.^2),2));
Normals=horzcat(Normals(:,1)./Distance,Normals(:,2)./Distance,Normals(:,3)./Distance);
f(:,4)=sqrt(sum((Normals-Normalsoriginal).^2,2));
f1=f(f(:,4)<1,1:3);
f2=f(f(:,4)>1,1:3);
f2(:,4)=f2(:,2);
f2(:,2)=[];
fnew=[fnew; f1; f2];
fnew=removerows(fnew,e(:,2));
[vnew]=project(vnew, fnew,voriginal,foriginal);
end
% [vnew, fnew]=cleanpatch(vnew, fnew);
% [vnew, fnew,temp] = edgecollaps( vnew, fnew, 0.1 );
end

function [vnew, fnew]=cleanpatch(V, F)
%remove duplicate vertices
 [vnew, indexm, indexn] =  unique(V, 'rows');
fnew = indexn(F);
%remove nonsens faces
numfaces = (1:size(fnew,1))';
e1=fnew(:,1)-fnew(:,2);
e2=fnew(:,1)-fnew(:,3);
e3=fnew(:,2)-fnew(:,3);
e1=[e1 numfaces];
e2=[e2 numfaces];
e3=[e3 numfaces];
e1=e1(e1(:,1)==0,2);
e2=e2(e2(:,1)==0,2);
e3=e3(e3(:,1)==0,2);
nonsensefaces=unique(vertcat(e1,e2,e3));
fnew=removerows(fnew,nonsensefaces);
% remove nonconnected vertices
numvertices = (1:size(vnew,1))';
connected=unique(reshape(fnew,3*size(fnew,1),1));
numvertices=removerows(numvertices,connected);
vtemp=vnew;
vtemp=removerows(vtemp,numvertices);
[lia,loc]=ismember(vnew,vtemp,'rows');
fnew = loc(fnew);
vnew=vtemp;
end
