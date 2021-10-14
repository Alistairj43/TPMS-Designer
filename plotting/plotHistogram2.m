function h = plotHistogram2(ax,FV,property1,property2,opts)
%% Function to calculate area weighted histograms
% Inputs:
% ax   - axis to plot on
% data - global data structure
% n - number of bin edges
%Outputs can be used in a plot/bar chart ex: plot(cax,y,yp);
% [x, yp, yc] - [Values of bin centres,...
%       Probability of a value falling within a bin,...
%       Cumulative probability distribution] 

pID = convertStringsToChars(extractBefore(property1+' ', ' '));
pID2 = convertStringsToChars(extractBefore(property2+' ', ' '));
if isfield(FV,pID)
    X = FV.(pID);
    Y = FV.(pID2);
elseif isfield(FV.Vproperty,pID)
    X = FV.Vproperty.(pID);
    Y = FV.Vproperty.(pID2);
elseif isfield(FV.Fproperty,pID)
    X = FV.Fproperty.(pID);
    Y = FV.Fproperty.(pID2);
end

if isfield(opts,'n')
    n = opts.n;
else
    n = 100;
end

h = histogram2(ax,X,Y,n,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
xlabel(ax,property1); ylabel(ax,property2);
c=colorbar(ax); c.Label.String = "Probability Density"; colormap(ax,hot);
end