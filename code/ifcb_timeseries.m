close all; clear all; clc;

cd features_t;

subdirout=dir(['*.csv']);
subdiroutcell=struct2cell(subdirout);
filelist=subdiroutcell(1,:);

for q = 1:(length(filelist));
f = csvread(char(filelist(q)),1,2);
v(1:size(f,1),q) = f(:,1);
l(1:size(f,1),q) = f(:,7);
q
end
%l = (6./pi.*v).^(1/3);
clearvars -EXCEPT v l;

%%

for i = 1:452;
    L = l(:,i);
    L = L(L>0);
    a(i) = plfit(L,'xmin',min(L));
    u(i) = plvar(L,'xmin',min(L),'reps',10);
    i
end

%%

A = sum(a./u.^2)./sum(u.^(-2))
U = sqrt(sum((A-a).^2.*(u.^(-2)))./sum((u.^(-2))))./sqrt(length(u))

cd ..