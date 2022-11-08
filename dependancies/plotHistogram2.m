function h = plotHistogram2(FV,property1,property2,opts,ax)
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
    FV SurfaceMesh = [];
    property1 string = [];
    property2 string = [];
    opts = [];
    ax = [];
end

if isempty(ax)
    figure; ax = gca;
end

pID = convertStringsToChars(extractBefore(property1+' ', ' '));
pID2 = convertStringsToChars(extractBefore(property2+' ', ' '));

if isfield(FV,pID)
    X = FV.(pID);
    Y = FV.(pID2);
elseif isfield(FV.Vproperty,pID)
    type = 'v';
    X = FV.Vproperty.(pID);
    Y = FV.Vproperty.(pID2);
elseif isfield(FV.Fproperty,pID)
    type = 'f';
    X = FV.Fproperty.(pID);
    Y = FV.Fproperty.(pID2);
else
    type = 'n';
    X = [];
    Y = [];
end

if isfield(opts,'n')
    n = max(ceil(opts.n+1),2);
else
    n = 100;
end

h = histogram2(ax,X,Y,n,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
xlabel(ax,property1); ylabel(ax,property2);
c=colorbar(ax); c.Label.String = "Probability Density"; colormap(ax,hot);
end