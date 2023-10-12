function [solid,U] = voxelateLattice(X,Y,Z,cellsize,nodes,struts,rstrut,rnode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n is the number of voxel along each axis
% address is the file location of wireframe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solid = zeros(size(X));      % initial grid with zeros
U = ones(size(X))*10000;
rstrut = rstrut./cellsize(1);
rnode = rnode./cellsize(1);

%% Get the voxel close the the strut witnin a certain distance
for i = 1:numel(X)          % for each voxel, deside if it is active
    % for each strut, get the distance to the voxel
    center = [X(i) Y(i) Z(i)]./cellsize-floor([X(i) Y(i) Z(i)]./cellsize); % voxel center position
    for j = 1:length(struts)
        start_n = nodes(struts(j,1),:);  % start node coordinate
        end_n = nodes(struts(j,2),:);    % end node coordinate
        % determine alpha and beta are acute angle
            distance1 = norm(cross(end_n - start_n,center - start_n))...
                /norm(end_n - start_n);
            distance2 = min(norm(center - start_n),norm(center - end_n));
        if (distance1-rstrut<=0||distance2-rnode<0)     % if distance less than radius, active it
            solid(i) = 1;
        end

        tempDist = min(distance1-rstrut,distance2-rnode);
        U(i) = min(U(i),tempDist);
    end
end
end
