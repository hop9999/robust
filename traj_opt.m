clear all;
clc;
close all;

dt = 1e-2;

k_i = 0.031;
mgl = 0.1852;
I = 0.0046;
k_q = 0.0012;
x_d = pi/2;

A = [0  1
    -mgl*cos(x_d)/I -k_q/I];
 
B = [0 
     k_i/I];
 
G = 0.1*[0
     mgl*sin(x_d)/I];

X_0 = 0.1*[2 -2
           1 1 ];

x0 = [-1 -1];

options = optimoptions('fmincon','Display','iter');
x = fmincon(@fun,x0,[],[],[],[],[],[],@constr,options)
dt = 1e-2;
x_lim = 1.1;
dx_lim = 1.1;
u_lim = 60;
x_f = 0.15;
dx_f = 0.05;

k_i = 0.03;
mgl = 0.183;
I = 0.0040;
k_q = 0.001;
x_d = pi/2;

A = [0  1 0  0
    -mgl*cos(x_d)/I -k_q/I 0 0
    0  0 0  1
    0  0 -mgl*cos(x_d)/I -k_q/I];

B = [0  0 
     k_i/I 0
    0  0 
     0 k_i/I];

G = 0.1*[0
     mgl*sin(x_d)/I
     0
     -mgl*sin(x_d)/I];
 


X_0 = 0.1*[2 -2 -2 2 
           -1 -1 1 1
           2 -2 -2 2 
           -1 -1 1 1];

K = [x(1) x(2) 0 0
     0 0 x(1) x(2)];

hold on
pgon = polyshape(X_0(1,:),X_0(2,:));
plot(pgon)
A_d = A*dt + [eye(2),zeros(2,2);zeros(2,2),eye(2)];
B_d = B*dt;
G_d = G*dt;

for i = 1:100
    X_0 = A_d*X_0 + G_d + B_d*K*X_0;
    


    if length(X_0) > 16
          X_n1 = interp1( 1:1:length(X_0),X_0(1,:),linspace(1,length(X_0),N_v));
         X_n2 = interp1( 1:1:length(X_0),X_0(2,:),linspace(1,length(X_0),N_v));
         X_n3 = interp1( 1:1:length(X_0),X_0(3,:),linspace(1,length(X_0),N_v));
         X_n4 = interp1( 1:1:length(X_0),X_0(4,:),linspace(1,length(X_0),N_v));
         X_0 = [X_n1;X_n2;X_n3;X_n4];
    end
    pgon = polyshape(X_0(1,:),X_0(2,:));
    plot(pgon)
    pgon = polyshape(X_0(3,:),X_0(4,:));
    plot(pgon)
end
plot(X_0(1,:),X_0(2,:),"*")
plot(X_0(3,:),X_0(4,:),"*")