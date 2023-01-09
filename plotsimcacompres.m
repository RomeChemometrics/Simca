function plotsimcacompres(smod, figmer, ncl)
eval(['x=smod.CMres{',num2str(ncl), '}.mean', figmer,';']);
eval(['sx=smod.CMres{',num2str(ncl), '}.std', figmer,';']);

M=1.1*max(x+sx); 

for i=1:4
    for j=1:3
       X{j}(i,:)=x(18*(i-1)+[2*j-1 2*j 2*j+5 2*j+6 2*j+11 2*j+12]);
       S{j}(i,:)=sx(18*(i-1)+[2*j-1 2*j 2*j+5 2*j+6 2*j+11 2*j+12]);
       
    end
end
X{4}=x(73:74); 
S{4}=sx(73:74); 
ngroups = 4;
nbars = 6;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
w=0.7*groupwidth/nbars; 




figure('units', 'normalized', 'position', [0.05 0.05 0.9 0.9])
subplot(2,2,1)
b1=bar(X{1});
b1(1).EdgeColor=[153 0 0]/255; b1(1).FaceColor=[153 0 0]/255;
b1(2).EdgeColor=[255 102 102]/255;b1(2).FaceColor=[255 102 102]/255;
b1(3).EdgeColor=[255 128 0]/255; b1(3).FaceColor=[255 128 0]/255; 
b1(4).EdgeColor=[255 178 102]/255;b1(4).FaceColor=[255 178 102]/255;
b1(5).EdgeColor=[204 204 0]/255; b1(5).FaceColor=[204 204 0]/255;
b1(6).EdgeColor=[255 255 102]/255;b1(6).FaceColor=[255 255 102]/255;

for i = 1:nbars
    xx = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    for j=1:ngroups
    RCerrorbar1(xx(j), X{1}(j,i)-S{1}(j,i), X{1}(j,i)+S{1}(j,i), w);
    end
end
set(gca, 'fontsize', 16, 'fontweight', 'bold', 'linewidth', 2)
set(gca, 'xticklabel', {'Percentile', 'F-distr (approx.)', 'F-distr (rig.)', '\chi^2'}, 'xticklabelrotation',325)
ylim([0 M])
ylabel(figmer)
title('Sym-SIMCA')

subplot(2,2,2)
b2=bar(X{2});
b2(1).EdgeColor=[153 0 0]/255; b2(1).FaceColor=[153 0 0]/255;
b2(2).EdgeColor=[255 102 102]/255;b2(2).FaceColor=[255 102 102]/255;
b2(3).EdgeColor=[255 128 0]/255; b2(3).FaceColor=[255 128 0]/255; 
b2(4).EdgeColor=[255 178 102]/255;b2(4).FaceColor=[255 178 102]/255;
b2(5).EdgeColor=[204 204 0]/255; b2(5).FaceColor=[204 204 0]/255;
b2(6).EdgeColor=[255 255 102]/255;b2(6).FaceColor=[255 255 102]/255;

for i = 1:nbars
    xx = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    for j=1:ngroups
    RCerrorbar1(xx(j), X{2}(j,i)-S{1}(j,i), X{2}(j,i)+S{1}(j,i), w);
    end
end
set(gca, 'fontsize', 16, 'fontweight', 'bold', 'linewidth', 2)
set(gca, 'xticklabel', {'Percentile', 'F-distr (approx.)', 'F-distr (rig.)', '\chi^2'}, 'xticklabelrotation',325)
ylim([0 M])
ylabel(figmer)
legend({'Perc.(rig)','Perc.(compl)', 'J.-M. (rig)','J.-M. (compl)', '\chi^2 (rig)', '\chi^2 (compl)'}, 'fontsize',13, 'fontweight', 'bold', 'position',[0.4675 0.8238 0.0810 0.1011])
title('Alt-SIMCA')

subplot(2,2,3)
b3=bar(X{3});
b3(1).EdgeColor=[153 0 0]/255; b3(1).FaceColor=[153 0 0]/255;
b3(2).EdgeColor=[255 102 102]/255;b3(2).FaceColor=[255 102 102]/255;
b3(3).EdgeColor=[255 128 0]/255; b3(3).FaceColor=[255 128 0]/255; 
b3(4).EdgeColor=[255 178 102]/255;b3(4).FaceColor=[255 178 102]/255;
b3(5).EdgeColor=[204 204 0]/255; b3(5).FaceColor=[204 204 0]/255;
b3(6).EdgeColor=[255 255 102]/255;b3(6).FaceColor=[255 255 102]/255;

for i = 1:nbars
    xx = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    for j=1:ngroups
    RCerrorbar1(xx(j), X{3}(j,i)-S{1}(j,i), X{3}(j,i)+S{1}(j,i), w);
    end
end
set(gca, 'fontsize', 16, 'fontweight', 'bold', 'linewidth', 2)
set(gca, 'xticklabel', {'Percentile', 'F-distr (approx.)', 'F-distr (rig.)', '\chi^2'}, 'xticklabelrotation',325)
ylim([0 M])
ylabel(figmer)
title('Combined Index')

subplot(2,2,4)
bar(1,X{4}(1), 'facecolor',[102 0 204]/255, 'edgecolor', [102 0 204]/255); 
hold on
bar(2,X{4}(2), 'facecolor',[255 51 153]/255, 'edgecolor', [255 51 153]/255); 
RCerrorbar1(1:2, X{4}-S{4}, X{4}+S{4}, 0.4);
set(gca, 'xtick', 1:2, 'xticklabel', {'Rigorous', 'Compliant'}, 'xticklabelrotation',325)
 
set(gca, 'fontsize', 16, 'fontweight', 'bold', 'linewidth', 2)
ylim([0 M])
ylabel(figmer)
title('Data Driven (Pomerantsev)')
