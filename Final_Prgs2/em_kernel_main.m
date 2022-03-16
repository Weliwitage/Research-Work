clear all
close all

load sample_y y;
load sample_e1 e1;
load sample_e2 e2;
load sample_e3 e3;

x=y;
clear y;
N=length(x);

P1=zeros(1,N);
P2=zeros(1,N);
P3=zeros(1,N);

f1=zeros(1,N);
f2=zeros(1,N);
f3=zeros(1,N);

MAXIT=10;

LL1=zeros(1,MAXIT);
LL2=zeros(1,MAXIT);
LL3=zeros(1,MAXIT);

PP1=zeros(N,MAXIT);
PP2=zeros(N,MAXIT);
PP3=zeros(N,MAXIT);

% Bandwidth 
h=1;

% Initial values 

L1=1.0/3;
L2=1.0/3;
L3=1.0/3;

sg1=1;
sg2=1;
sg3=1;

mu1=-3;
mu2=0;
mu3=3;

for i=1:N
  f1(i)=(1/sg1)*(1/(sqrt(2*pi)))*exp(-0.5*((x(i)-mu1)/sg1)^2);
  f2(i)=(1/sg2)*(1/(sqrt(2*pi)))*exp(-0.5*((x(i)-mu2)/sg2)^2);
  f3(i)=(1/sg3)*(1/(sqrt(2*pi)))*exp(-0.5*((x(i)-mu3)/sg3)^2);
end
for i=1:N
    P1(i)=L1*f1(i)/(L1*f1(i)+L2*f2(i)+L3*f3(i));
    P2(i)=L2*f2(i)/(L1*f1(i)+L2*f2(i)+L3*f3(i));
    P3(i)=L3*f3(i)/(L1*f1(i)+L2*f2(i)+L3*f3(i)); 
end
for q=1:10
L1=sum(P1)/N;
L2=sum(P2)/N;
L3=sum(P3)/N;

sumf1=zeros(1,N);
sumf2=zeros(1,N);
sumf3=zeros(1,N);

for j=1:N
   for k=1:N 
     sumf1(i)=sumf1(i)+P1(k)*kernel((x(k)-x(i))/h);
     sumf2(i)=sumf1(i)+P2(k)*kernel((x(k)-x(i))/h);
     sumf3(i)=sumf1(i)+P3(k)*kernel((x(k)-x(i))/h);
   end
end
f1(i)=sumf1(i)/(N*h*L1);
f2(i)=sumf2(i)/(N*h*L2);
f3(i)=sumf3(i)/(N*h*L3);

LL1(q)=L1;
LL2(q)=L2;
LL3(q)=L3;

PP1(:,q)=P1(:);
PP2(:,q)=P2(:);
PP3(:,q)=P3(:);
end

D=[PP1(:,10) PP2(:,10) PP3(:,10)];

C=zeros(N,1);
I=zeros(N,1);
set1=0;
set2=0;
set3=0;

for i=1:N
    [C(i),I(i)]=max(D(i,:));
    if (I(i)==1) 
        set1=set1+1;
    elseif(I(i)==2)
        set2=set2+1;
    else 
        set3=set3+1;
    end
end

d1=zeros(set1,1);
d2=zeros(set2,1);
d3=zeros(set3,1);
k1=0;
k2=0;
k3=0;

for i=1:N
    if (I(i)==1) 
       k1=k1+1;d1(k1)=x(i); 
    elseif (I(i)==2)
       k2=k2+1;d2(k2)=x(i);
    else
       k3=k3+1; d3(k3)=x(i);
    end
end

figure(10)
plot(e1,'g.');hold on
plot(e2,'r.');
plot(e3,'b.'); hold on

figure(11)
plot(d1,'bo','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerSize',10);
            hold on
plot(d2,'ro','LineWidth',2,...
                'MarkerEdgeColor','r',...
                'MarkerSize',10);
plot(d3,'go','LineWidth',2,...
                'MarkerEdgeColor','b',...
                'MarkerSize',10);
title('After clustering via EM algorithm with nonparametric');

hold off



