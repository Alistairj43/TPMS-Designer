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

if isfield(data.FV,pID) % FV General Property
    x_in = data.FV.(pID);
elseif isfield(data.FV.Fproperty,pID) % Face Property
    x_in = data.FV.Fproperty.(pID);
    try
        a_in = data.FV.Fproperty.area;
    catch
        a_in = ones(size(x_in));
    end
elseif isfield(data.FV.Vproperty,pID) % Vertex Property
    x_in = data.FV.Vproperty.(pID);
    a_in = ones(size(x_in));
elseif isfield(data.F.property,pID) % Field Property
    
    x_in = data.F.property.(pID);
    [l,w,h] = size(x_in);
    x_in = reshape(x_in,[l*w*h,1,1]); % Flatten
    a_in = ones(size(x_in));
else
    x_in = [0 1];
    a_in = ones(size(x_in));
    pName = "Data does not exist";
end

if isfield(opts,'n')
    n = max(ceil(opts.n+1),2);
else
    n = 100;
end


temp = prctile(x_in,[0.2 99.8]);
try
    binEdges = linspace(temp(1),temp(2),n);
catch
    binEdges = linspace(0,1,n);
end

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