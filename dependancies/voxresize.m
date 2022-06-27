function out = voxresize(in,sz)

[l,w,h] = size(in);

if l==1 && w==1 && h==1
    out = in; % No reshape nessecary
else
    if l==1
        in(2,:,:) = in(1,:,:);
    end
    if w==1
        in(:,2,:) = in(:,1,:);
    end
    if h==1
        in(:,:,2) = in(:,:,1);
    end

    [y, x, z] = ndgrid(linspace(1,size(in,1),sz(1)),...
        linspace(1,size(in,2),sz(2)),...
        linspace(1,size(in,3),sz(3)));
    out=interp3(in,x,y,z); % Return the rehspaped array
end
