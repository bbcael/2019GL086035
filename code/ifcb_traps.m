close all; clear all; clc;

cd features_1and2;

subdirout=dir(['*150.csv']);
subdiroutcell=struct2cell(subdirout);
filelist=subdiroutcell(1,:);

for q = 1:(length(filelist));
f = csvread(char(filelist(q)),1,2);
v(1:size(f,1),q) = f(:,1);
l(1:size(f,1),q) = f(:,7);
q
end
l = (6./pi.*v).^(1/3);
clearvars -EXCEPT v l;

%L = l(l>0);
%a = plfit(L,'xmin',min(L));
%u = plvar(L,'xmin',min(L),'reps',100);

%%

for i = 1:size(l,2);
    L = l(:,i);
    L = L(L>0);
    a(i) = plfit(L,'xmin',3.667);
    u(i) = plvar(L,'xmin',3.667,'reps',100);
    i
end

%%

A = sum(a./u.^2)./sum(u.^(-2))
U = sqrt(sum((A-a).^2.*(u.^(-2)))./sum((u.^(-2))))./sqrt(length(u))

cd ..