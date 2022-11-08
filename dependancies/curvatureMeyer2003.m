function [GC, MC, k1, k2]= curvatureMeyer2003(x,y,z,tri3d)
% calculation of Gaussian (GC) and mean curvatures (MC) of a discrete surface. 
% The coordinates of the points are given by x, y, z.
% "tri" is a triangulation table which gives the vertex ID of each triangle. 
% The IDs is consistent with the order of data points in x,y,z. 
% It returns GC and MC at each point.
% Without loss of any generality, in notation, it was supposed that the
% points of each triangle are in counter-clockwise order from 1 to 3.
% By starting from vertex 1, the vectors of edges are also in
% counter-clockwise from v1 to v3.
% 
% The method is based on this paper:
% Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003). 
% Discrete differential-geometry operators for triangulated 2-manifolds. 
% In Visualization and mathematics III (pp. 35-57). Springer Berlin Heidelberg.
% Available at : http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.24.3427&rep=rep1&type=pdf
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: The curvatures data at the oundaries of domain is not reliable. It
% artificially gives zero instead of non-sense data.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Update version 1.1:
% The above paper, gives a vector for MC and the absoulute value of mean
% curvature is half of the norm of this vector. Therfore, in this way the sign is
% not given. Here, I used a dot product betwen the MC vector and the normal
% vector at each point calculated based on weighted averaging of the
% triangle normal vectors given by MATLAB. I'm not sure what is the convention of
% MATLAB in determining the direction for the normal vectors. But it seems
% that they are all consistent toward one side of the surface. So the
% calculated signed MC shows the change of sign in data, but it can not
% gaurantee that, for example, the positive MC is for a locally convex
% region (if The GC is positive).
% 
% Written by Alireza Dastan   V.1.1        21/01/2017
% email: ardastan@gmail.com

bndry_edge=freeBoundary(tri3d);
if isempty(bndry_edge)
    bndry_edge = 0;
end
f_normal = faceNormal(tri3d);
f_center = incenter(tri3d);
tri = tri3d.ConnectivityList;
%%% area
for i=1:length(tri(:,1))
    v1=[x(tri(i,2))-x(tri(i,1)),y(tri(i,2))-y(tri(i,1)),z(tri(i,2))-z(tri(i,1))];
    v2=[x(tri(i,3))-x(tri(i,1)),y(tri(i,3))-y(tri(i,1)),z(tri(i,3))-z(tri(i,1))];
    area_tri(i,1)=0.5*norm(cross(v1,v2));
end
%%%% angles and edges of each triangle
for i=1:length(tri(:,1))
     
    p1=tri(i,1);
    p2=tri(i,2);
    p3=tri(i,3);
    
    v1(i,:)=[x(p2)-x(p1),y(p2)-y(p1),z(p2)-z(p1)];
    v2(i,:)=[x(p3)-x(p2),y(p3)-y(p2),z(p3)-z(p2)];
    v3(i,:)=[x(p1)-x(p3),y(p1)-y(p3),z(p1)-z(p3)];
    
    l_edg(i,1)=norm(v1(i,:));
    l_edg(i,2)=norm(v2(i,:));
    l_edg(i,3)=norm(v3(i,:));
    
    ang_tri(i,1)=acos(dot(v1(i,:)/l_edg(i,1),-v3(i,:)/l_edg(i,3)));
    ang_tri(i,2)=acos(dot(-v1(i,:)/l_edg(i,1),v2(i,:)/l_edg(i,2)));
    ang_tri(i,3)=pi-(ang_tri(i,1)+ang_tri(i,2));
end
a_mixed=zeros(1,length(x));
alf=zeros(1,length(x));
GC=zeros(length(x),1);
MC=zeros(length(x),1);
for i=1:length(x)
    mc_vec=[0,0,0];
    n_vec=[0,0,0];
    
    if find(bndry_edge(:,1)==i)
    else
    clear neib_tri
    neib_tri=vertexAttachments(tri3d,i);
    
    for j=1:length(neib_tri{1})
        neib=neib_tri{1}(j);
        
        %%%% sum of angles around point i ===> GC
        for k=1:3
            if tri(neib,k)==i
              alf(i)=alf(i)+ ang_tri(neib,k);
              break;
            end
        end
        
        %%%%% mean curvature operator
        if     k==1
            mc_vec=mc_vec+(v1(neib,:)/tan(ang_tri(neib,3))-v3(neib,:)/tan(ang_tri(neib,2)));
        elseif k==2
            mc_vec=mc_vec+(v2(neib,:)/tan(ang_tri(neib,1))-v1(neib,:)/tan(ang_tri(neib,3)));
        elseif k==3
            mc_vec=mc_vec+(v3(neib,:)/tan(ang_tri(neib,2))-v2(neib,:)/tan(ang_tri(neib,1)));
        end
        
        
        %%% A_mixed calculation
        if(ang_tri(neib,k)>=pi/2)
            a_mixed(i)=a_mixed(i)+area_tri(neib)/2;
        else
            if (any(ang_tri(neib,:)>=pi/2))
                a_mixed(i)=a_mixed(i)+area_tri(neib)/4;
            else
                sum=0;
                for m=1:3
                    if m~=k
                        ll=m+1;
                        if ll==4       %% p1==>l2   ,p2==>l3   ,p3==>l1    
                            ll=1;
                        end
                        sum=sum+(l_edg(neib,ll)^2/tan(ang_tri(neib,m)));
                    end
                end
                a_mixed(i)=a_mixed(i)+sum/8;
            end
        end
        
        %%%% normal vector at each vertex    
        %%%% weighted average of normal vecotors of neighbour triangles
        wi=1/norm([f_center(neib,1)-x(i),f_center(neib,2)-y(i),f_center(neib,3)-z(i)]);
        n_vec=n_vec+wi*f_normal(neib,:);
        
    end
       
    GC(i)=(2*pi()-alf(i))/a_mixed(i);
    
    mc_vec=0.25*mc_vec/a_mixed(i);
    n_vec=n_vec/norm(n_vec);
    %%%% sign of MC
    if dot(mc_vec,n_vec) <0
        MC(i)=-norm(mc_vec);
    else
        MC(i)=norm(mc_vec);
    end
    
    end
end

k1 = real(MC+sqrt(MC.^2-GC));
k2 = real(MC-sqrt(MC.^2-GC));
end
