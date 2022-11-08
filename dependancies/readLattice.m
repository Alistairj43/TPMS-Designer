%%Import information of strut
function [nodes,struts] = readLattice(filepathname)
fid = fopen(filepathname,'r');
k = 1;
j = 1;
tline = fgetl(fid);
while ischar(tline)
    if (tline(1) == 'G')
        x = str2double(tline(17:24));
        y = str2double(tline(25:32));
        z = str2double(tline(33:40));
        nodes(k,1:3) = [x,y,z];
        k = k + 1;
    end
    if (tline(1) == 'S')
        Snode = str2double(tline(17:24));
        Enode = str2double(tline(25:32));
        struts(j,1:2) = [Snode,Enode];
        j = j + 1;
    end
    tline = fgetl(fid);
end
fclose(fid);
end