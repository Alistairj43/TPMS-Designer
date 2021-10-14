% input the homogenized elasticity tensor
function s = plotTensor(ax,CH)
% transform it to 3*3*3*3 tensor
tensor = generate(CH);
% find the E1 in 360 degree
% x = 0:pi/180:2*pi;
[a,e] = meshgrid(0:0.05*pi:2*pi, -pi/2:0.025*pi:pi/2);
E1 = zeros(size(a));
for i = 1:size(a,1)
    for j = 1:size(a,2)
        % build the transformation matrix
        trans_z = [cos(a(i,j))     -sin(a(i,j))      0;
                   sin(a(i,j))      cos(a(i,j))      0;
                   0                0                1];
        trans_y = [cos(e(i,j))      0           sin(e(i,j));
                   0                1           0;
                  -sin(e(i,j))      0           cos(e(i,j))];
        % calculate the new tensor
        N_tensor = transform(tensor,trans_y*trans_z);
        % transform the tensor to 6*6
        N_CH = ToMatrix(N_tensor);
        % calculate the E1
        E = modulus(N_CH);
        E1(i,j) = E(1);
    end
end
[x,y,z] = sph2cart(a,e,E1);
c = sqrt(x.^2+y.^2+z.^2);
s = surf(ax,x,y,z,c,'FaceColor','interp','EdgeColor','none');
c=colorbar(ax); c.Label.String = 'Elastic Modulus'; colormap(ax,"jet");
axis(ax,'manual','vis3d','equal','tight');
xlabel(ax,"X"); ylabel(ax,"Y"); zlabel(ax,"Z");
ax.XTick = []; ax.YTick = []; ax.ZTick = [];
ax.BoxStyle = 'full'; 
rotate3d(ax,'on');
end
function [E] = modulus(CH)
S = inv(CH);
E = zeros(6,1);
E(1) = 1/S(1,1);
E(2) = 1/S(2,2);
E(3) = 1/S(3,3);
E(4) = 1/S(4,4);
E(5) = 1/S(5,5);
E(6) = 1/S(6,6);
end
function C = generate(CH)
C = zeros(3,3,3,3);
for i = 1:6
    for j = 1:6
        [a,b] = change(i);
        [c,d] = change(j);
        C(a,b,c,d) = CH(i,j);
    end
end
for i = 1:3
    if (i == 3)
        j = 1;
    else
        j = i+1;
    end
    for m = 1:3
        if (m == 3)
            n = 1;
        else
            n = m+1;
        end
        C(j,i,n,m) = C(i,j,m,n);
        C(j,i,m,n) = C(i,j,m,n);
        C(i,j,n,m) = C(i,j,m,n);
        C(j,i,m,m) = C(i,j,m,m);
        C(m,m,j,i) = C(m,m,i,j);
    end
end
end
% change the index 4 5 6 to 23 31 12
function [a,b] = change(w)
if (w < 4)
    a = w;
    b = w;
else
    if (w == 4)
        a = 2;
        b = 3;
    else
        if (w == 5)
            a = 3;
            b = 1;
        else
            if (w==6)
                a = 1;
                b = 2;
            end
        end
    end
end
end
function CH = ToMatrix(C)
CH = zeros(6,6);
for i = 1:6
    for j = 1:6
        [a,b] = change(i);
        [c,d] = change(j);
        CH(i,j) = C(a,b,c,d);     
    end
end
end

function otr = transform(itr,tmx)
%
% FUNCTION
% otr = transform(itr,tmx)
%
% DESCRIPTION
% transform 3D-tensor (Euclidean or Cartesion tensor) of any order (>0) to another coordinate system
%
% PARAMETERS
% otr = output tensor, after transformation; has the same dimensions as the input tensor
% itr = input tensor, before transformation; should be a 3-element vector, a 3x3 matrix, or a 3x3x3x... multidimensional array, each dimension containing 3 elements
% tmx = transformation matrix, 3x3 matrix that contains the direction cosines between the old and the new coordinate system
%
ne = numel(itr);                % number of tensor elements
nd = ndims(itr);                % number of tensor dimensions, i.e. order of tensor
if (ne==3), nd = 1; end         % order of tensor is 1 in case of a 3x1 or 1x3 vector
otr = itr;                      % create output tensor
otr(:) = 0;                     % fill output tensor with zeros; this way a symbolic tensor remains symbolic
cne = cumprod(3*ones(nd,1))/3;  % vector with cumulative number of elements for each dimension (divided by three)
for oe = 1:ne                                 % loop over all output elements
   ioe = mod(floor((oe-1)./cne),3)+1;          % calculate indices of current output tensor element
   for ie = 1:ne                             % loop over all input elements
      pmx = 1;                                 % initialise product of transformation matrices
      iie = mod(floor((ie-1)./cne),3)+1;       % calculate indices of current input tensor element
      for id = 1:nd                          % loop over all dimensions
         pmx = pmx * tmx( ioe(id), iie(id) );  % create product of transformation matrices
      end
      otr(oe) = otr(oe) + pmx * itr(ie);       % add product of transformation matrices and input tensor element to output tensor element
   end
end
% Transform matrix about Z axis
% for x = 0:pi/20:pi/2
%     trans = [cos(x)      cos(x-pi/2) 0;
%              cos(x+pi/2) cos(x)      0;
%              0           0           1];
%     N_CH = transform(C,trans);
% end
end