clear all
close all

N = 600;
k = 3;
MAXIT=100;

% Definition of Matrices 

cp=zeros(MAXIT,3);
cm=zeros(MAXIT,3);
cs=zeros(MAXIT,3);

y  = zeros(1,N);
p  = zeros(1,3);
f  = zeros(N,k);
w  = zeros(N,k);
sg = zeros(1,k);
mu = zeros(1,k);
mu1= zeros(N,k);

% % Generation of a mixture of three normally distributed data
% load data_e1 e1;
% load data_e2 e2;
% load data_e3 e3;

% e1 = normrnd(-4,2.0,[1,N/3]);
% e2 = normrnd(0,2.0,[1,(N/3)]);
% e3 = normrnd(4,2.0,[1,(N/3)]);

T1=[0.5 0.3 0.2 -5 0 2 2 0.5 1];
save t1_data T1;

%   for q=1:N/3
%      y(3*q-2)=e1(q);
%      y(3*q-1)=e2(q);
%      y(3*q)  =e3(q);
%   end
%   
% save data_e1 e1;
% save data_e2 e2;
% save data_e3 e3;
% 

%y=[e1,e2,e3];
%Second method of generating data set 

y=[];N=600;
for i=1:N
    ra=rand(1,1);
    add=2*randn-5;
    if ra <0.2 
        add=randn+2;
    elseif ra >0.7
        add=0.5*randn;
    end
y=[y add];
end
save sample_y y;
clear y;
% 
% save sample_e1 e1;
% save sample_e2 e2;
% save sample_e3 e3;

load sample_y y;

% Initial Values
p(1,1:k) = 1/k;
sg(1,1:k)= 1;
mu(1,1)= -4;
mu(1,2)= 1;
mu(1,3)= 1;

% Em Algorithm 
for r=1:MAXIT
    for j=1:k
        for i=1:N
            f(i,j)=(1/sg(j))*(1/(sqrt(2*pi)))*exp(-0.5*((y(i)-mu(j))/sg(j))^2);
            w(i,j)=(f(i,j)*p(j))/(f(i,1)*p(1)+f(i,2)*p(2)+f(i,3)*p(3));
        end    
    end
    p  = sum(w,1)/N;
    mu =y*w./sum(w,1);
    mu1 = repcolumn(mu',N)';
    for j=1:k
        temp1=0;
        temp2=0;
        for i=1:N
            temp1=temp1+(w(i,j)*(y(i)-mu1(i,j))^2);
            temp2=temp2+w(i,j);
        end
        sg(j)=sqrt(temp1/temp2);
    end
    
    % For EM iteration steps    
cp(r,1:3)=p(1:3);
cm(r,1:3)=mu(1:3);    
cs(r,1:3)=sg(1:3);   
end

T2=[cp(MAXIT,1:3) cm(MAXIT,1:3) cs(MAXIT,1:3)];
save t2_data T2;

Z=zeros(N,k);
for i=1:N
  for j=1:k
     Z(i,j)=cp(MAXIT,j)*(1.0/sqrt(2*pi*cs(MAXIT,j)))*...
            exp(-(y(i)-cm(MAXIT,j))^2/(2*cs(MAXIT,j)^2));
  end
end
C=zeros(N,1);
I=zeros(N,1);
set1=0;
set2=0;
set3=0;
for i=1:N
    [C(i),I(i)]=max(Z(i,:));
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
       k1=k1+1;d1(k1)=y(i); 
    elseif (I(i)==2)
       k2=k2+1;d2(k2)=y(i);
    else
       k3=k3+1; d3(k3)=y(i);
    end
end

save data_d1 d1;
save data_d2 d2;
save data_d3 d3;

% Generation of points for evaluation of the density function

N=length(y);
x = min(y):(max(y)-min(y))/(N/2):max(y);
x=x';
z=zeros(length(x),1);

% Nonparametric Density estimation 
h =1.3;
[f] = main_regression(x,z,y,h);

% Evaluation of explicit function using the estimations of mixture model

for k=1:length(x)
   z(k)=p(1)*(1/(sg(1)*sqrt(2*pi)))*exp((-(x(k)-mu(1))^2/(2*sg(1)*sg(1))))+...
        p(2)*(1/(sg(2)*sqrt(2*pi)))*exp((-(x(k)-mu(2))^2/(2*sg(2)*sg(2))))+...
        p(3)*(1/(sg(3)*sqrt(2*pi)))*exp((-(x(k)-mu(3))^2/(2*sg(3)*sg(3))));  
end

% % Representations by graphs
figure(1)
  plot(x,f,'r-','LineWidth',2);
  hold on
  plot(x,z,'g-','LineWidth',2);
  legend('kernel','Gaussian Mixture');
  
figure(2)
hold on
plot(cp(:,1),'g-','LineWidth',2);
plot(cp(:,2),'r-','LineWidth',2);
plot(cp(:,3),'b-','LineWidth',2);
hold off
legend('p1','p2','p3');
title('Convergence of P values');
xlabel('Iteration Step');
ylabel('p-value');
grid on

figure(3)
hold on
plot(cm(:,1),'g-','LineWidth',2);
plot(cm(:,2),'r-','LineWidth',2);
plot(cm(:,3),'b-','LineWidth',2);
hold off
legend('mu1','mu2','mu3');
title('Convergence of means');
xlabel('Iteration Step');
ylabel('mu-value');
grid on

figure(4)
hold on
plot(cs(:,1),'g-','LineWidth',2);
plot(cs(:,2),'r-','LineWidth',2);
plot(cs(:,3),'b-','LineWidth',2);
hold off
legend('sg1','sg2','sg3');
title('Convergence of variance');
xlabel('Iteration Step');
ylabel('sigma-value');
grid on

figure(5)
plot(y,'k*');
title('Mixture of 3 normally distributed data set');

figure(6)
plot(y,'k*');
hold on 
for i=1:N
    if (I(i)==1)
    plot(i,y(i),'go','LineWidth',2,...
                'MarkerEdgeColor','b',...
                'MarkerSize',10);
    elseif (I(i)==2)
    plot(i,y(i),'go','LineWidth',2,...
                'MarkerEdgeColor','r',...
                'MarkerSize',10);   
    else 
    plot(i,y(i),'go','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerSize',10);
    end
end
title('Model based clustering via EM algorithm');
hold off



%%% Clustering via Nonparametric Method with EM

load sample_y y;
% load sample_e1 e1;
% load sample_e2 e2;
% load sample_e3 e3;

x=y;
clear y;
N=length(x);

P1=zeros(1,N);
P2=zeros(1,N);
P3=zeros(1,N);

f1=zeros(1,N);
f2=zeros(1,N);
f3=zeros(1,N);

MAXIT=40;

LL1=zeros(1,MAXIT);
LL2=zeros(1,MAXIT);
LL3=zeros(1,MAXIT);

PP1=zeros(N,MAXIT);
PP2=zeros(N,MAXIT);
PP3=zeros(N,MAXIT);

% Bandwidth 
h=((4.0/3)^(1.0/5))*(N)^(-1.0/5);

% Initial values 

L1=1.0/3;
L2=1.0/3;
L3=1.0/3;

sg1=1;
sg2=1;
sg3=1;

mu1=-4;
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
for q=1:MAXIT
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

D=[PP1(:,MAXIT) PP2(:,MAXIT) PP3(:,MAXIT)];

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

s1=zeros(set1,1);
s2=zeros(set2,1);
s3=zeros(set3,1);

k1=0;
k2=0;
k3=0;

for i=1:N
    if (I(i)==1) 
       k1=k1+1;s1(k1)=x(i); 
    elseif (I(i)==2)
       k2=k2+1;s2(k2)=x(i);
    else
       k3=k3+1; s3(k3)=x(i);
    end
end
save data_s1 s1;
save data_s2 s2;
save data_s3 s3;


figure(8)
plot(x,'k*');
hold on 
for i=1:N
    if (I(i)==1)
    plot(i,x(i),'go','LineWidth',2,...
                'MarkerEdgeColor','r',...
                'MarkerSize',10);
    elseif (I(i)==2)
    plot(i,x(i),'go','LineWidth',2,...
                'MarkerEdgeColor','b',...
                'MarkerSize',10);   
    else 
    plot(i,x(i),'go','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerSize',10);
    end
end
title('After clustering via EM algorithm with nonparametric');
hold off

[length(d1),length(s1);length(d2),length(s2);length(d3),length(s3)]
[mean(d1),mean(s1);mean(d2),mean(s2);mean(d3),mean(s3)]
[var(d1),var(s1);var(d2),var(s2);var(d3),var(s3)]



