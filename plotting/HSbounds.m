plots = ["elastic", "poisson", "zenerRatio", "shear", "totalStiffness"];
plotlabels = ["Elastic Modulus", "Poisson Ratio", "Zener Anisotropy Ratio", "Shear Modulus", "Total Stiffness"];
colors = [0 0 0; 1 0 0; 0 0 1;0 0 0; 1 0 0; 0 0 1];
sym = '...xxx';

%% HS Bounds
Gb = 0.3759; Ga = 0;
Eb = 1; Ea = 0;
vb = 0.33; va = 0;
Kb = Eb/(3-6*vb); Ka = Ea/(3-6*va);
Vb = 0:0.01:1;
Va = 1-Vb;
K_upper = Kb+Va./((1/(Ka-Kb))+(3*Vb/(3*Kb+4*Gb)));
G_upper = Gb+Va./((1/(Ga-Gb))+(6*Vb*(Kb+2*Gb)/(5*Gb*(3*Kb+4*Gb))));
E1 = 2*G_upper.*(1+vb);
E2 = 3*K_upper.*(1-2*vb);
HS.bulk = Kb*4*Gb*(Vb)./(4*Gb+3*Kb.*Va);
HS.shear = Gb*(9*Kb+8*Gb)*Vb./(20*Gb+15*Kb-6*Vb*(Kb+2*Gb));
HS.elastic = Eb*9*(HS.shear).*(HS.bulk)./((3*HS.bulk)+(HS.shear));
HS.poisson = (3*HS.bulk-2*HS.shear)./(6*HS.bulk+2*HS.shear);
HS.totalStiffness = HS.elastic+2.*HS.shear.*(1-HS.poisson);
HS.zenerRatio = HS.elastic./(2*HS.shear.*(1+HS.poisson));

%% Plotting
for i=1:length(plots)
    f = figure;
    ax = gca;
    L = gscatter(ax,DOEout.volumeFraction,DOEout.(plots(i)),strcat(DOEout.equation, DOEout.type),colors,sym);
    hold on
    plot(ax,Vb,HS.(plots(i)),'DisplayName','HS-Upper Bound')
    xlabel('Volume Fraction');
    ylabel(plotlabels(i));
    grid on;
end