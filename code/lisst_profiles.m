clear all; close all; clc;
% profile: slope, volume-mean diameter, integrated volume

load('lisstKM1910_L1700759_sphere.mat')

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

V = sum(corr_vd,2);

scatter(smooth(depth(is1:im),V(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),V(im:is2),20),depth(im:is2),10,'filled')

z1 = depth(is1:im);
v1 = V(is1:im);
z2 = depth(im:is2);
v2 = V(im:is2);

load('lisstKM1910_L1681056_sphere.mat')

corr_vd = corr_vd(1:2450,:);
depth = depth(1:2450);

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

V = sum(corr_vd,2);

scatter(smooth(depth(is1:im),V(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),V(im:is2),20),depth(im:is2),10,'filled')

z3 = depth(is1:im);
v3 = V(is1:im);
z4 = depth(im:is2);
v4 = V(im:is2);

load('lisstKM1910_L1681056_sphere.mat')

corr_vd = corr_vd(2451:end,:);
depth = depth(2451:end);

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

V = sum(corr_vd,2);

scatter(smooth(depth(is1:im),V(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),V(im:is2),20),depth(im:is2),10,'filled')

z5 = depth(is1:im);
v5 = V(is1:im);
z6 = depth(im:is2);
v6 = V(im:is2);

plot(300:301,300:301)
plot(.9206.*exp((20:240).*-0.01431),20:240)
plot(.4527.*exp((20:240).*-0.01381),20:240)
plot(1.387.*exp((20:240).*-0.02904),20:240)
plot(.3733.*exp((20:240).*-0.01391),20:240)
plot(.6439.*exp((20:240).*-0.01761),20:240)
plot(.3641.*exp((20:240).*-0.01297),20:240)


lgd = legend('down 1','up 1','down 2a','up 2a','down 2b','up 2b','location','southeast');
set(lgd,'interpreter','latex','fontsize',20)
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('integrated particle volume [$\mu$L/L]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
box on
axis([0 .82 0 250])

%%

clear all; close all; clc;
% profile: slope, volume-mean diameter, integrated volume

load('lisstKM1910_L1700759_sphere.mat')

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

M = sum(dias.*corr_vd,2)./sum(corr_vd,2);

scatter(smooth(depth(is1:im),M(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),M(im:is2),20),depth(im:is2),10,'filled')

z1 = depth(is1:im);
m1 = M(is1:im);
z2 = depth(im:is2);
m2 = M(im:is2);

load('lisstKM1910_L1681056_sphere.mat')

corr_vd = corr_vd(1:2450,:);
depth = depth(1:2450);

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

M = sum(dias.*corr_vd,2)./sum(corr_vd,2);

scatter(smooth(depth(is1:im),M(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),M(im:is2),20),depth(im:is2),10,'filled')

z3 = depth(is1:im);
m3 = M(is1:im);
z4 = depth(im:is2);
m4 = M(im:is2);

load('lisstKM1910_L1681056_sphere.mat')

corr_vd = corr_vd(2451:end,:);
depth = depth(2451:end);

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

M = sum(dias.*corr_vd,2)./sum(corr_vd,2);

scatter(smooth(depth(is1:im),M(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),M(im:is2),20),depth(im:is2),10,'filled')

z5 = depth(is1:im);
m5 = M(is1:im);
z6 = depth(im:is2);
m6 = M(im:is2);

%plot(300:301,300:301)
%plot(.9206.*exp((20:240).*-0.01431),20:240)
%plot(.4527.*exp((20:240).*-0.01381),20:240)
%plot(1.387.*exp((20:240).*-0.02904),20:240)
%plot(.3733.*exp((20:240).*-0.01391),20:240)
%plot(.6439.*exp((20:240).*-0.01761),20:240)
%plot(.3641.*exp((20:240).*-0.01297),20:240)


lgd = legend('down 1','up 1','down 2a','up 2a','down 2b','up 2b','location','southeast');
set(lgd,'interpreter','latex','fontsize',20)
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('volume-mean particle diameter [$\mu$m]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
box on
axis([0 210 0 250])

%%

clear all; close all; clc;
% profile: slope, volume-mean diameter, integrated volume

load('lisstKM1910_L1700759_sphere.mat')

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

n = 6./pi.*corr_vd./dias.^3;

for i = 1755:(size(corr_vd,1)-24);
    mn = nanmean(n(i-25:i+24,:),1);
    wn = 1./nanvar(n(i-25:i+24,:),1);
    modelFun = @(b,x) b(1).*x.^(-b(2));
    start = [1,5];
    wnlm = fitnlm(dias(9:27),mn(9:27),modelFun,start,'Weight',wn(9:27));
    xi(i) = wnlm.Coefficients{2,1};
    xi_pm(i) = wnlm.Coefficients{2,2};
    i
end

%%

%xi(1756:1846) = NaN;
xi(xi>10) = NaN;

scatter(smooth(depth(is1:im),xi(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),xi(im:is2),20),depth(im:is2),10,'filled')

z1 = depth(is1:im);
x1 = xi(is1:im);
z2 = depth(im:is2);
x2 = xi(im:is2);

%%

load('lisstKM1910_L1681056_sphere.mat')

corr_vd = corr_vd(1:2450,:);
depth = depth(1:2450);

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

n = 6./pi.*corr_vd./dias.^3;

for i = 1283:(size(corr_vd,1)-24);
    mn = nanmean(n(i-25:i+24,:),1);
    wn = 1./nanvar(n(i-25:i+24,:),1);
    modelFun = @(b,x) b(1).*x.^(-b(2));
    start = [.1,4.5];
    wnlm = fitnlm(dias(9:27),mn(9:27),modelFun,start,'Weight',wn(9:27));
    xi(i) = wnlm.Coefficients{2,1};
    xi_pm(i) = wnlm.Coefficients{2,2};
    i
end

xi(2421:2450) = NaN;
xi(xi>10) = NaN;

%%

scatter(smooth(depth(is1:im),xi(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),xi(im:is2),20),depth(im:is2),10,'filled')

z3 = depth(is1:im);
x3 = xi(is1:im);
z4 = depth(im:is2);
x4 = xi(im:is2);

%%

load('lisstKM1910_L1681056_sphere.mat')

corr_vd = corr_vd(2451:end,:);
depth = depth(2451:end);

[m,im] = max(depth);
is = find(depth<20);
is1 = max(is(is<im));
is2 = min(is(is>im));
clear m;

n = 6./pi.*corr_vd./dias.^3;

for i = 1442:(size(corr_vd,1)-24);
    mn = nanmean(n(i-25:i+24,:),1);
    wn = 1./nanvar(n(i-25:i+24,:),1);
    modelFun = @(b,x) b(1).*x.^(-b(2));
    start = [.1,4.5];
    wnlm = fitnlm(dias(9:27),mn(9:27),modelFun,start,'Weight',wn(9:27));
    xi(i) = wnlm.Coefficients{2,1};
    xi_pm(i) = wnlm.Coefficients{2,2};
    i
end

xi(2238:2276) = NaN;
xi(xi>10) = NaN;

%%

scatter(smooth(depth(is1:im),xi(is1:im),20),depth(is1:im),10,'filled')
set(gca,'ydir','reverse')
hold on
scatter(smooth(depth(im:is2),xi(im:is2),20),depth(im:is2),10,'filled')

z5 = depth(is1:im);
x5 = xi(is1:im);
z6 = depth(im:is2);
x6 = xi(im:is2);

%%

plot(300:301,300:301)
plot(.9206.*exp((20:240).*-0.01431),20:240)
plot(.4527.*exp((20:240).*-0.01381),20:240)
plot(1.387.*exp((20:240).*-0.02904),20:240)
plot(.3733.*exp((20:240).*-0.01391),20:240)
plot(.6439.*exp((20:240).*-0.01761),20:240)
plot(.3641.*exp((20:240).*-0.01297),20:240)


lgd = legend('down 1','up 1','down 2a','up 2a','down 2b','up 2b','location','southeast');
set(lgd,'interpreter','latex','fontsize',20)
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('integrated particle volume [$\mu$L/L]','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
box on
axis([0 .82 0 250])

%%

scatter(x1,z1,10,'filled')
set(gca,'ydir','reverse')
hold on
box on
scatter(x2,z2,10,'filled')
scatter(x3,z3,10,'filled')
scatter(x4,z4,10,'filled')
scatter(x5,z5,10,'filled')
scatter(x6,z6,10,'filled')

lgd = legend('down 1','up 1','down 2a','up 2a','down 2b','up 2b','location','northwest');
set(lgd,'interpreter','latex','fontsize',20)
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('power-law exponent','Interpreter','latex','FontSize',20)
ylabel('depth [m]','Interpreter','latex','FontSize',20)
axis([1 4 0 250])

