close all
clear all


load sample_y y;

% original data with means 
%load data_e1 e1;
%load data_e2 e2;
%load data_e3 e3;

% model based clustering 
load data_d1 d1; d1=d1';
load data_d2 d2; d2=d2';
load data_d3 d3; d3=d3';

% nonparametric clustering 
load data_s1 s1; s1=s1';
load data_s2 s2; s2=s2';
load data_s3 s3; s3=s3';


% Load parameters 

load t1_data T1; T1=T1';
load t2_data T2; T2=T2';

% Generation of distributions from given data set with correct parameters
x=-8:0.1:8;
fe1=zeros(1,length(x));
fe2=zeros(1,length(x));
fe3=zeros(1,length(x));

for k=1:length(x)
  fe1(k)=(1.0/(T1(7)*sqrt(2*pi)))*exp((-(x(k)-T1(4))^2/(2*T1(7)*T1(7))));
  fe2(k)=(1.0/(T1(8)*sqrt(2*pi)))*exp((-(x(k)-T1(5))^2/(2*T1(8)*T1(8))));
  fe3(k)=(1.0/(T1(9)*sqrt(2*pi)))*exp((-(x(k)-T1(6))^2/(2*T1(9)*T1(9))));
end

figure(1)
plot(x,fe1,'r-','LineWidth',2);
hold on 
plot(x,fe2,'r-','LineWidth',2);
plot(x,fe3,'r-','LineWidth',2);


% Generation of distributions from gaussian mixture 
x=-8:0.1:8;
fd1=zeros(1,length(x));
fd2=zeros(1,length(x));
fd3=zeros(1,length(x));

for k=1:length(x)
  fd1(k)=(1/(T2(7)*sqrt(2*pi)))*exp((-(x(k)-T2(4))^2/(2*T2(7)*T2(7))));
  fd2(k)=(1/(T2(8)*sqrt(2*pi)))*exp((-(x(k)-T2(5))^2/(2*T2(8)*T2(8))));
  fd3(k)=(1/(T2(9)*sqrt(2*pi)))*exp((-(x(k)-T2(6))^2/(2*T2(9)*T2(9))));
end
pause(3);

figure(1)
plot(x,fd1,'b-','LineWidth',2);
hold on 
plot(x,fd2,'b-','LineWidth',2);
plot(x,fd3,'b-','LineWidth',2);



% Generation of points for evaluation of the density function
% Nonparametric Density estimations

h=1.5;

N1=length(s1);
x1 = min(s1):(max(s1)-min(s1))/(N1/2):max(s1);
x1=x1';
x1=-8:0.1:8;
z1=zeros(length(x1),1);
[fn1] = main_regression(x1,z1,s1,h);

N2=length(s2);
x2 = min(s2):(max(s2)-min(s2))/(N2/2):max(s2);
x2=x2';
x2=-8:0.1:8;
z2=zeros(length(x2),1);
[fn2] = main_regression(x2,z2,s2,h);

N3=length(s3);
x3 = min(s3):(max(s3)-min(s3))/(N3/2):max(s3);
x3=x3';
x3=-8:0.1:8;
z3=zeros(length(x3),1);
[fn3] = main_regression(x3,z3,s3,h);

pause(3);

figure(1)
plot(x1,fn1,'g-','LineWidth',2);
hold on 
plot(x2,fn2,'g-','LineWidth',2);
plot(x3,fn3,'g-','LineWidth',2);
hold off
grid on
%legend('exact','parametric','nonparametric');
