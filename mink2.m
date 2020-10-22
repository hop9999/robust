clear all;
clc;
close all;

dt = 1e-2;

k_i = 0.03;
mgl = 0.183;

I = 0.0040;
k_q = 0.001;

x_d = pi/2;

A_min = [0                1
         -mgl*cos(x_d)/I  -0.9*k_q/I];

A_max = [0                1
         -mgl*cos(x_d)/I  -10.9*k_q/I];

B_min = [0 
         k_i/I];
 
G_min = [0
         -0.1*mgl*sin(x_d)/I];
     
G_max = [0
         0.1*mgl*sin(x_d)/I];
     
A = [A_max zeros(2,2)
     zeros(2,2) A_min];

A = [A zeros(4,4)
     zeros(4,4) A];
 
B = [B_min zeros(2,1)
     zeros(2,1) B_min];
 
B = [B zeros(4,2)
     zeros(4,2) B]
 
G = [G_min
     G_min
     G_max
     G_max];
 
M = [ones(1,4) zeros(1,4) zeros(1,4) zeros(1,4)
     zeros(1,4) ones(1,4) zeros(1,4) zeros(1,4)
     zeros(1,4) zeros(1,4) ones(1,4) zeros(1,4)
     zeros(1,4) zeros(1,4) zeros(1,4) ones(1,4)];

x_0 = [2 -2 -2 2 
       -1 -1 1 1]
   
X_0 = [x_0
       x_0
       x_0
       x_0]*M;
   
K = [-10 -1]  
K = [zeros(1,0) K zeros(1,6)
     zeros(1,2) K zeros(1,4)
     zeros(1,4) K zeros(1,2)
     zeros(1,6) K zeros(1,0)];

hold on
pgon = polyshape(X_0(1,:),X_0(2,:));
plot(pgon)
A_d = A*dt + [eye(8)];
          
B_d = B*dt;
G_d = G*dt;

for i = 1:100

    X_0 = A_d*X_0 + G_d + B_d*K*X_0;
    length(X_0)
    if length(X_0) > 16
         X_n1 = interp1( 1:1:length(X_0),X_0(1,:),linspace(1,length(X_0),N_v));
         X_n2 = interp1( 1:1:length(X_0),X_0(2,:),linspace(1,length(X_0),N_v));
         X_n3 = interp1( 1:1:length(X_0),X_0(3,:),linspace(1,length(X_0),N_v));
         X_n4 = interp1( 1:1:length(X_0),X_0(4,:),linspace(1,length(X_0),N_v));
         X_n5 = interp1( 1:1:length(X_0),X_0(5,:),linspace(1,length(X_0),N_v));
         X_n6 = interp1( 1:1:length(X_0),X_0(6,:),linspace(1,length(X_0),N_v));
         X_n7 = interp1( 1:1:length(X_0),X_0(7,:),linspace(1,length(X_0),N_v));
         X_n8 = interp1( 1:1:length(X_0),X_0(8,:),linspace(1,length(X_0),N_v));
         X_0 = [X_n1;X_n2;X_n3;X_n4;X_n5;X_n6;X_n7;X_n8];
    end
    pgon = polyshape(X_0(1,:),X_0(2,:));
    plot(pgon)
    pgon = polyshape(X_0(3,:),X_0(4,:));
    plot(pgon)
    pgon = polyshape(X_0(5,:),X_0(6,:));
    plot(pgon)
    pgon = polyshape(X_0(7,:),X_0(8,:));
    plot(pgon)
end
X_0
plot(X_0(1,:),X_0(2,:),"*")
plot(X_0(3,:),X_0(4,:),"*")
plot(X_0(5,:),X_0(6,:),"*")
plot(X_0(7,:),X_0(8,:),"*")