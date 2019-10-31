close all; clear all; clc;

% first section: create mat files with aggregate properties

subdirout=dir(['*2*300*.mat']);
subdiroutcell=struct2cell(subdirout);
filelist=subdiroutcell(1,:);

for q = 1:(length(filelist)./3);
    clearvars -EXCEPT q filelist D M N S V;
    close all

load(filelist{(3*(q-1)+1)});
disp(filelist{(3*(q-1)+1)});
l = 1e-6.*dias;
v1 = 1e-6.*corr_vd(1:20,:);
mn1 = sum(l.*v1,2)./sum(v1,2);
load(filelist{(3*(q-1)+2)});
disp(filelist{(3*(q-1)+2)});
v2 = 1e-6.*corr_vd(1:20,:);
mn2 = sum(l.*v2,2)./sum(v2,2);
load(filelist{(3*(q-1)+3)});
disp(filelist{(3*(q-1)+3)});
v3 = 1e-6.*corr_vd(1:20,:);
mn3 = sum(l.*v3,2)./sum(v3,2);
v = 1e6.*[sum(v1,2); sum(v2,2); sum(v3,2)];
m = [mn1; mn2; mn3];
clearvars -EXCEPT v1 v2 v3 v l m D M N S V q filelist;

%% change to number density

n1 = 6./pi.*v1./l.^3;
n2 = 6./pi.*v2./l.^3;
n3 = 6./pi.*v3./l.^3;
d = [n1; n2; n3];

for i = 1:size(d,1); % slope the lazy way
[a s(i) aa] = regression(log10(l),log10(d(i,:)));
end
clear a aa;

n = sum(d,2);

clear i n1 n2 n3 v1 v2 v3;

D(:,:,q) = d;
M(:,q) = m;
N(:,q) = n;
S(:,q) = s;
V(:,q) = v;

end

clearvars -EXCEPT D l M N S V;

%save 2_300.mat;
%plot(l,d)
%set(gca,'xscale','log','yscale','log')
%plot(diasP,NdP,'k','linewidth',2);
%set(gca,'yscale','log','xscale','log')
%hold on
%plot(diasP,NdP+Nd_errP,'--','color',[.59 .07 .39],'linewidth',2)
%axis([1 250 1e-6 0.3])
%set(gca,'TickLabelInterpreter','latex','FontSize',20)
%xlabel('Equivalent Spherical Diameter [$\mu$m]','Interpreter','latex','FontSize',20)
%ylabel('Probability Density [$\mu$m$^{-1}$]','Interpreter','latex','FontSize',20)
%title('First Deployment, 300m, Trap L, Replicate 3','Interpreter','latex','FontSize',20)
%box on

%%

%{
clear all; close all; clc;

plot(117.*exp(-(65:310)./345.6),65:310,'-.k','linewidth',2)
hold on

load 2_75.mat;
S = -S;
M = 1e6.*M;
V(42,8) = NaN; % first deployment
[q1,ind] = sort(nanmean(M,1));
qq1 = nanstd(M,1);
errorbar(q1,linspace(70,80,11),qq1(ind),'horizontal','+k')
hold on
errorbar(nanmean(M(:)),75,nanstd(M(:)),'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

m75 = nanmean(M(:)); s75 = nanstd(M(:));

load 2_150.mat;
S = -S;
M = 1e6.*M;
%V([1 2],1) = NaN; % second deployment
%V(3,5) = NaN; % second deployment
[q2,ind] = sort(nanmean(M,1));
qq2 = nanstd(M,1);
errorbar(q2,linspace(145,155,11),qq2(ind),'horizontal','+k')
errorbar(nanmean(M(:)),150,nanstd(M(:)),'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)
axis([0 20 69 170])
set(gca,'ydir','reverse')

m150 = nanmean(M(:)); s150 = nanstd(M(:));

load 2_300.mat;
S = -S;
M = 1e6.*M;
V(1,6) = NaN; % first deployment
%V(1,5) = NaN; % second deployment
[q3,ind] = sort(nanmean(M,1));
qq3 = nanstd(M,1);
errorbar(q3,linspace(295,305,11),qq3(ind),'horizontal','+k')
errorbar(nanmean(M(:)),300,nanstd(M(:)),'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)
axis([10 125 65 310])
set(gca,'ydir','reverse')

m300 = nanmean(M(:)); s300 = nanstd(M(:));

set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('volume-mean particle diameter [$\mu$m]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('second deployment','Interpreter','latex','FontSize',20)
box on
%}

zz = [75 150 300]; mm = [m75 m150 m300]; ss = [s75 s150 s300].^(-2); vv = mm.^3; ssv = [s75 s150 s300].^(-6);

%text(20,85,'$\ell_V \approx 110$m','interpreter','latex','fontsize',20)
%}

%%

clear all; close all; clc;

%plot(8.327.*exp(-(65:310).*.007342),65:310,'-.','linewidth',2,'color',[.75 .75 .75])
plot(15.38.*exp(-(65:310).*0.009762),65:310,'-.k','linewidth',2,'color',[.75 .75 .75])
hold on

load 2_75.mat;
%V(V>30) = NaN; % first deployment, one outlier
[q1,ind] = sort(nanmean(V,1));
qq1 = nanstd(V,1);
errorbar(q1,linspace(70,80,11),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,75,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

m75 = nanmean(V(:)); s75 = nanstd(V(:));

load 2_150.mat;
V(V>30) = NaN; %second deployment, 3 outliers
[q2,ind] = sort(nanmean(V,1));
qq2 = nanstd(V,1);
errorbar(q2,linspace(145,155,11),qq2(ind),'horizontal','+k')
wm2 = sum(q2./(qq2.^2))./sum((qq2.^(-2)))
ws2 = sqrt(sum((q2-wm2).^2.*(qq2.^(-2)))./sum((qq2.^(-2))))./sqrt(length(q2))
errorbar(wm2,150,ws2,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

m150 = nanmean(V(:)); s150 = nanstd(V(:));

load 2_300.mat;
V(V>10) = NaN; % first+second deployment, one outlier eat
[q3,ind] = sort(nanmean(V,1));
qq3 = nanstd(V,1);
errorbar(q3,linspace(295,305,11),qq3(ind),'horizontal','+k')
wm3 = sum(q3./(qq3.^2))./sum((qq3.^(-2)))
ws3 = sqrt(sum((q3-wm3).^2.*(qq3.^(-2)))./sum((qq3.^(-2))))./sqrt(length(q3))
errorbar(wm3,300,ws3,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)
set(gca,'ydir','reverse')

set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('integrated particle volume [$\mu$L/L]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('second deployment','Interpreter','latex','FontSize',20)
box on
axis([0 15 65 310])

%}

zz = [75 150 300]; ww = [wm wm2 wm3]; ss = [ws ws2 ws3].^(-2); 

%text(20,85,'$\ell_V \approx 110$m','interpreter','latex','fontsize',20)

%%

clear all; close all; clc;
load 2_75.mat;

D(V>30) = NaN;

for i = 1:11;  
ds = squeeze(D(:,7:27,i));
dm = nanmean(ds,1);
dw = 1./var(ds,1);
sl = l(7:27);
modelFun = @(b,x) b(1).*x.^(-b(2));
start = [sum(dm),2];
wnlm = fitnlm(sl,dm,modelFun,start,'Weight',dw);
xi_75_1(i) = wnlm.Coefficients{2,1};
xi_pm_75_1(i) = wnlm.Coefficients{2,2};
i
end

load 2_150.mat;

D(V>30) = NaN;

for i = 1:11;  
ds = squeeze(D(:,7:27,i));
dm = nanmean(ds,1);
dw = 1./nanvar(ds,1);
sl = l(7:27);
modelFun = @(b,x) b(1).*x.^(-b(2));
start = [sum(dm),2];
wnlm = fitnlm(sl,dm,modelFun,start,'Weight',dw);
xi_150_1(i) = wnlm.Coefficients{2,1};
xi_pm_150_1(i) = wnlm.Coefficients{2,2};
i
end

load 2_300.mat;

D(V>10) = NaN;

for i = 1:11;  
ds = squeeze(D(:,7:27,i));
dm = nanmean(ds,1);
dw = 1./var(ds,1);
sl = l(7:27);
modelFun = @(b,x) b(1).*x.^(-b(2));
start = [sum(dm),2];
wnlm = fitnlm(sl,dm,modelFun,start,'Weight',dw);
xi_300_1(i) = wnlm.Coefficients{2,1};
xi_pm_300_1(i) = wnlm.Coefficients{2,2};
i
end

clearvars -EXCEPT xi*;

%%

clear all; close all; clc;
load xi_1.mat;
load xi_2.mat;

%plot(14.65.*exp(-(65:310)./103.3),65:310,'-.k','linewidth',2)
%hold on

[q1,ind] = sort(xi_75_2);
qq1 = xi_pm_75_2(ind);;
errorbar(q1,linspace(70,80,11),qq1,'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,75,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

[q1,ind] = sort(xi_150_2);
qq1 = xi_pm_150_2(ind);;
errorbar(q1,linspace(145,155,11),qq1,'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,150,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

set(gca,'ydir','reverse')

[q1,ind] = sort(xi_300_2);
qq1 = xi_pm_300_2(ind);;
errorbar(q1,linspace(295,305,11),qq1,'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,300,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('power-law exponent','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('first deployment','Interpreter','latex','FontSize',20)
box on
axis([1.75 2.55 65 310])

%}

%zz = [75 150 300]; mm = [m75 m150 m300]; ss = [s75 s150 s300].^(-2);

%text(20,85,'$\ell_V \approx 110$m','interpreter','latex','fontsize',20)


% NB:
% n(1K75_2, 2C75_1, 2K150_2 = 19, 20th is avg of first 19. 2H150_3 is average of 2H150_1 & 2H150_2)

%%


clear all; close all; clc;

plot(10.47.*exp(-(65:310)./117.3),65:310,'-.k','linewidth',2)
%plot(14.65.*exp(-(65:310)./103.3),65:310,'-.k','linewidth',2)
hold on

load 1_75.mat;
V(V>30) = NaN; % first deployment, one outlier
[q1,ind] = sort(nanmean(V,1));
qq1 = nanstd(V,1);
errorbar(q1,linspace(70,80,12),qq1(ind),'horizontal','+k')
hold on
errorbar(nanmean(V(:)),75,nanstd(V(:)),'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

m75 = nanmean(V(:)); s75 = nanstd(V(:));

load 1_150.mat;
%V(V>30) = NaN; %second deployment, 3 outliers
[q2,ind] = sort(nanmean(V,1));
qq2 = nanstd(V,1);
errorbar(q2,linspace(145,155,12),qq2(ind),'horizontal','+k')
errorbar(nanmean(V(:)),150,nanstd(V(:)),'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)
%axis([0 20 69 170])
set(gca,'ydir','reverse')

m150 = nanmean(V(:)); s150 = nanstd(V(:));

load 1_300.mat;
V(V>10) = NaN; % first+second deployment, one outlier eat
[q3,ind] = sort(nanmean(V,1));
qq3 = nanstd(V,1);
errorbar(q3,linspace(295,305,12),qq3(ind),'horizontal','+k')
errorbar(nanmean(V(:)),300,nanstd(V(:)),'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)
set(gca,'ydir','reverse')

m300 = nanmean(V(:)); s300 = nanstd(V(:));

set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('integrated particle volume [$\mu$L/L]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('first deployment','Interpreter','latex','FontSize',20)
box on
axis([0 17 65 310])

%}

zz = [75 150 300]; mm = [m75 m150 m300]; ss = [s75 s150 s300].^(-2);

%text(20,85,'$\ell_V \approx 110$m','interpreter','latex','fontsize',20)

%%

load 1_75.mat;

L = repmat(l,60,1,12);
m = squeeze(sum(L(:,7:27,:).*D(:,7:27,:),2)./sum(D(:,7:27,:),2))*1e6;
[q1,ind] = sort(nanmean(m,1));
qq1 = nanstd(m,1);
errorbar(q1,linspace(70,80,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,75,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)


load 1_150.mat;

L = repmat(l,60,1,12);
m = squeeze(sum(L(:,7:27,:).*D(:,7:27,:),2)./sum(D(:,7:27,:),2))*1e6;
[q1,ind] = sort(nanmean(m,1));
qq1 = nanstd(m,1);
errorbar(q1,linspace(145,155,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,150,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

load 1_300.mat;

L = repmat(l,60,1,12);
m = squeeze(sum(L(:,7:27,:).*D(:,7:27,:),2)./sum(D(:,7:27,:),2))*1e6;
[q1,ind] = sort(nanmean(m,1));
qq1 = nanstd(m,1);
errorbar(q1,linspace(295,305,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,300,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

set(gca,'TickLabelInterpreter','latex','FontSize',16,'ydir','reverse')
xlabel('mean particle diameter [m]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('second deployment','Interpreter','latex','FontSize',20)
box on
axis([0 1e-5 65 310])

%%

clear all; close all; clc;
load 1_75.mat;

V = 1e-6.*V;
d = 1e6.*6./pi.*(V./N).^(1/3);
[q1,ind] = sort(nanmean(d,1));
qq1 = nanstd(d,1);
errorbar(q1,linspace(70,80,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,75,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

load 1_150.mat;

V = 1e-6.*V;
d = 1e6.*6./pi.*(V./N).^(1/3);
[q1,ind] = sort(nanmean(d,1));
qq1 = nanstd(d,1);
errorbar(q1,linspace(145,155,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,150,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

load 1_300.mat;

V(V>10) = NaN; % first+second deployment, one outlier eat
V = 1e-6.*V;
d = 1e6.*6./pi.*(V./N).^(1/3);
[q1,ind] = sort(nanmean(d,1));
qq1 = nanstd(d,1);
errorbar(q1,linspace(295,305,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,300,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)


set(gca,'TickLabelInterpreter','latex','FontSize',16,'ydir','reverse')
xlabel('volume-mean particle diameter [$\mu$m]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('first deployment','Interpreter','latex','FontSize',20)
box on
axis([0 30 65 310])


%%


L = repmat(l,60,1,11);
m = squeeze(sum(L(:,5:27,:).*D(:,5:27,:),2)./sum(D(:,5:27,:),2));
[q1,ind] = sort(nanmean(m,1));
qq1 = nanstd(m,1);
errorbar(q1,linspace(70,80,11),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,75,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)


load 2_150.mat;

L = repmat(l,60,1,11);
m = squeeze(sum(L(:,5:27,:).*D(:,5:27,:),2)./sum(D(:,5:27,:),2));
[q1,ind] = sort(nanmean(m,1));
qq1 = nanstd(m,1);
errorbar(q1,linspace(145,155,11),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,150,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

load 2_300.mat;

L = repmat(l,60,1,11);
m = squeeze(sum(L(:,5:27,:).*D(:,5:27,:),2)./sum(D(:,5:27,:),2));
[q1,ind] = sort(nanmean(m,1));
qq1 = nanstd(m,1);
errorbar(q1,linspace(295,305,11),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,300,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

set(gca,'TickLabelInterpreter','latex','FontSize',16,'ydir','reverse')
xlabel('mean particle diameter [m]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('second deployment','Interpreter','latex','FontSize',20)
box on
axis([0 1e-5 65 310])

%%

clear all; close all; clc;
load 1_75.mat;

[q1,ind] = sort(nanmean(N,1));
qq1 = nanstd(N,1);
errorbar(q1,linspace(70,80,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,75,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

load 1_150.mat;

N(N>3e10) = NaN;
[q1,ind] = sort(nanmean(N,1));
qq1 = nanstd(N,1);
errorbar(q1,linspace(145,155,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,150,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

load 1_300.mat;

[q1,ind] = sort(nanmean(N,1));
qq1 = nanstd(N,1);
errorbar(q1,linspace(295,305,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,300,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)


set(gca,'TickLabelInterpreter','latex','FontSize',16,'ydir','reverse')
xlabel('particle density [m$^{-3}$]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('first deployment','Interpreter','latex','FontSize',20)
box on
axis([0 Inf 65 310])

%%

clear all; close all; clc;

plot(.0001174.*exp(-(65:310).*.002966),65:310,'-.','linewidth',2,'color',[.75 .75 .75])
plot(.0001116.*exp(-(65:310).*.002966),65:310,'-.','linewidth',2,'color',[.75 .75 .75])
hold on

load 1_75.mat;
%V(V>30) = NaN; % first deployment, one outlier
[q1,ind] = sort(nanmean(M,1));
qq1 = nanstd(M,1);
errorbar(q1,linspace(70,80,12),qq1(ind),'horizontal','+k')
hold on
wm = sum(q1./(qq1.^2))./sum((qq1.^(-2)))
ws = sqrt(sum((q1-wm).^2.*(qq1.^(-2)))./sum((qq1.^(-2))))./sqrt(length(q1))
errorbar(wm,75,ws,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

m75 = nanmean(V(:)); s75 = nanstd(V(:));

load 1_150.mat;
V(V>30) = NaN; %second deployment, 3 outliers
[q2,ind] = sort(nanmean(M,1));
qq2 = nanstd(M,1);
errorbar(q2,linspace(145,155,12),qq2(ind),'horizontal','+k')
wm2 = sum(q2./(qq2.^2))./sum((qq2.^(-2)))
ws2 = sqrt(sum((q2-wm2).^2.*(qq2.^(-2)))./sum((qq2.^(-2))))./sqrt(length(q2))
errorbar(wm2,150,ws2,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)

m150 = nanmean(V(:)); s150 = nanstd(V(:));

load 1_300.mat;
V(V>10) = NaN; % first+second deployment, one outlier eat
[q3,ind] = sort(nanmean(M,1));
qq3 = nanstd(M,1);
errorbar(q3,linspace(295,305,12),qq3(ind),'horizontal','+k')
wm3 = sum(q3./(qq3.^2))./sum((qq3.^(-2)))
ws3 = sqrt(sum((q3-wm3).^2.*(qq3.^(-2)))./sum((qq3.^(-2))))./sqrt(length(q3))
errorbar(wm3,300,ws3,'horizontal','ok','color',[.59 .07 .39],'markersize',15,'markerfacecolor',[.59 .07 .39],'linewidth',3)
set(gca,'ydir','reverse')

set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('volume-mean particle diameter [m]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
title('first deployment','Interpreter','latex','FontSize',20)
box on
axis([0 .00011 65 310])

%}

zz = [75 150 300]; ww = [wm wm2 wm3]; ss = [ws ws2 ws3].^(-2); 

%text(20,85,'$\ell_V \approx 110$m','interpreter','latex','fontsize',20)