%%Import information of strut
function h = plotLattice(nodes,struts,ax)
% Function to plot lattice onto axis 'ax'
%Defaults
arguments
    nodes double;
    struts double = [];
    ax = [];
end

if isempty(ax)
    figure; ax = gca;
end

if ~isempty(struts)
    S = nodes(struts(:,1),:);
    E = nodes(struts(:,2),:);
    X = [S(:,1) E(:,1)]';
    Y = [S(:,2) E(:,2)]';
    Z = [S(:,3) E(:,3)]';
    h = plot3(ax,X,Y,Z,'-o','MarkerFaceColor','black');
elseif ~isempty(nodes)
    scatter3(ax,nodes(:,1),nodes(:,2), nodes(:,3))
end

axis(ax,'manual','vis3d','equal','tight');


end