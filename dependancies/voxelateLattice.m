function [voxel] = voxelateLattice(X,Y,Z,bounds,address,rstrut,rnode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n is the number of voxel along each axis
% address is the file location of wireframe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voxel = zeros(size(X));      % initial grid with zeros

%% Get the voxel close the the strut witnin a certain distance
[node,strut] = ReadStrut(address); % get the information of strut
node = node.*(bounds(2,:)-bounds(1,:))+bounds(1,:);
for i = 1:numel(X)          % for each voxel, deside if it is active
    % for each strut, get the distance to the voxel
    center = [X(i) Y(i) Z(i)];        % voxel center position
    for j = 1:length(strut)
        start_n = node(strut(j,1),:);  % start node coordinate
        end_n = node(strut(j,2),:);    % end node coordinate
        % determine alpha and beta are acute angle
            distance1 = norm(cross(end_n - start_n,center - start_n))...
                /norm(end_n - start_n);
            distance2 = min(norm(center - start_n),norm(center - end_n));
        if (distance1<=rstrut||distance2<rnode)     % if distance less than radius, active it
            voxel(i) = 1;
            continue;
        end
    end
end
end
%%Import information of strut
function [nodelist,strutlist] = ReadStrut(address)
fid = fopen(address,'r');
k = 1;
j = 1;
tline = fgetl(fid);
while ischar(tline)
    if (tline(1) == 'G')
        x = str2double(tline(17:24));
        y = str2double(tline(25:32));
        z = str2double(tline(33:40));
        nodelist(k,1:3) = [x,y,z];
        k = k + 1;
    end
    if (tline(1) == 'S')
        Snode = str2double(tline(17:24));
        Enode = str2double(tline(25:32));
        strutlist(j,1:2) = [Snode,Enode];
        j = j + 1;
    end
    tline = fgetl(fid);
end
fclose(fid);
end
