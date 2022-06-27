function h = plotHistogram(ax,FV,pName,opts)
%% Function to calculate area weighted histograms
% Inputs:
% ax   - axis to plot on
% data - global data structure
% n - number of bin edges
%Outputs can be used in a plot/bar chart ex: plot(cax,y,yp);
% [x, yp, yc] - [Values of bin centres,...
%       Probability of a value falling within a bin,...
%       Cumulative probability distribution] 

pID = convertStringsToChars(extractBefore(pName+' ', ' '));

if isfield(FV,pID)
    x_in = FV.(pID);
elseif isfield(FV.Fproperty,pID)
    x_in = FV.Fproperty.(pID);
    try
        a_in = FV.Fproperty.area;
    catch
        a_in = ones(size(x_in));
    end
elseif isfield(FV.Vproperty,pID)
    x_in = FV.Vproperty.(pID);
    a_in = ones(size(x_in));
else
    x_in = [0 1];
    pName = "Data does not exist";
    a_in = ones(size(x_in));
end


if isfield(opts,'n')
    n = opts.n;
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