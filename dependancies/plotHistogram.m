function h = plotHistogram(data,pName,opts,ax)
%% Function to calculate area weighted histograms
% Inputs:
% ax   - axis to plot on
% data - global data structure
% n - number of bin edges
%Outputs can be used in a plot/bar chart ex: plot(cax,y,yp);
% [x, yp, yc] - [Values of bin centres,...
%       Probability of a value falling within a bin,...
%       Cumulative probability distribution]
arguments
    data = [];
    pName string = [];
    opts = [];
    ax = [];
end

if isempty(ax)
    figure; ax = gca; view(3);
end

pID = convertStringsToChars(extractBefore(pName+' ', ' '));

if isfield(opts,'n')
    n = max(ceil(opts.n+1),2);
else
    n = 100;
end

try
    if isa(data,'UnitCell')
        if isfield(data.FV.Fproperty,pID)
            % Face property
            x_in = data.FV.Fproperty.(pID);
            try a_in = data.FV.Fproperty.area;
            catch; a_in = ones(size(x_in)); % Face areas not found
            end
        elseif isfield(data.FV.Vproperty,pID)
            % Vertex Property
            x_in = data.FV.Vproperty.(pID);
            a_in = ones(size(x_in));
        else
            % F property
            x_in = data.F.property.(pID);
            [l,w,h] = size(x_in);
            x_in = reshape(x_in,[l*w*h,1,1]); % Flatten
            a_in = ones(size(x_in));
        end
    elseif isa(data,'surfaceMesh')
        if isfield(data.Fproperty,pID)
            % Face property
            x_in = data.Fproperty.(pID);
            try a_in = data.Fproperty.area;
            catch; a_in = ones(size(x_in)); % Face areas not found
            end
        else
            % Vertex Property
            x_in = data.Vproperty.(pID);
            a_in = ones(size(x_in));
        end
    else % Couldn't determine type, assume raw data
        x_in = data;
        [l,w,h] = size(x_in);
        x_in = reshape(x_in,[l*w*h,1,1]); % Flatten
        a_in = ones(size_x_in);
    end
catch
    % Failed to validly find dataset
    x_in = [0 1];
    a_in = ones(size(x_in));
    pName = "Property data not found for: " + pName;
    n = 2;
end

temp = prctile(x_in,[0.2 99.8]);

if temp(1)==temp(2) % For datasets with no range
    temp(2) = temp(2)+0.5;
    temp(1) = temp(1)-0.5;
    n = 2;
    pName = "Data may not exist for: " + pName;
end

binEdges = linspace(temp(1),temp(2),n);

x = (binEdges(2:end)+binEdges(1:end-1))/2;
totalA = sum(a_in,"all","omitnan");
yp = zeros(size(x)); yc = zeros(size(x));
for i = 1:length(x)
    inds = x_in>=binEdges(i)&x_in<binEdges(i+1);
    yp(i) = sum(a_in(inds==1),"all","omitnan")/totalA;
    yc(i) = sum(yp,"all","omitnan");
end
yp = yp.*n;

% Do the plotting
cla(ax,'reset');
h(1) = bar(ax,x,yp); grid(ax,"on"); ylabel(ax,"Probability Density");
hold(ax,'on'); yyaxis(ax,'right'); h(2) = plot(ax,x,yc,'-r');
xlabel(ax,pName); ylabel(ax,"Cumulative Distribution");
hold(ax,"off");
end