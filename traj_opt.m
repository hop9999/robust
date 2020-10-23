clear all;
clc;
close all;

dt = 1e-2;

%model params
k_i = 0.03;
mgl = 0.183;
I = 0.0040;
k_q = 0.001;

x_d = pi/2;

%%%%%%%%%%%%%%%%%%%
%A nominal
A_c = [0                1
       -mgl*cos(x_d)/I  -5*k_q/I];
%A band
A_d1 = [0  0
        0  -5*k_q/I];
   
A_d2 = [0                    0
        -0.1*mgl*cos(x_d)/I  0];
    
%%%%%%%%%%%%%%%%%%%%
%B nominal
B_c = [0 
       k_i/I];
%B band
B_d = [0 
       0.1*k_i/I];
   
%%%%%%%%%%%%%%%%%%%
%G nominal
G_c = [0
       0];

%G band
G_d = [0
       0.1*mgl*sin(x_d)/I];
   
%%%%%%%%%%%%%%%%%%%%
% Matrix manipulations
n_comb = 16;

%Create matrix with all possible combinations of band params
A = blkdiag(kron(eye(4),A_c + A_d1),...
            kron(eye(4),A_c - A_d1),...
            kron(eye(4),A_c + A_d2),...
            kron(eye(4),A_c - A_d2));

B = blkdiag(kron(eye(2),B_c + B_d),...
            kron(eye(2),B_c - B_d),...
            kron(eye(2),B_c + B_d),...
            kron(eye(2),B_c - B_d),...
            kron(eye(2),B_c + B_d),...
            kron(eye(2),B_c - B_d),...
            kron(eye(2),B_c + B_d),...
            kron(eye(2),B_c - B_d));

G = [G_c - G_d
     G_c + G_d
     G_c - G_d
     G_c + G_d
     G_c - G_d
     G_c + G_d
     G_c - G_d
     G_c + G_d
     G_c - G_d
     G_c + G_d
     G_c - G_d
     G_c + G_d
     G_c - G_d
     G_c + G_d
     G_c - G_d
     G_c + G_d];

%init gain
K_0 = [-1 -1];

%init coordinates
N_init_x = 4;
x_0 = [0.1 -0.1 -0.1 0.1 
       -0.1 -0.1 0.1 0.1];
   
X_0 = kron(ones(n_comb,1),x_0);
   
%Dynamic discretization
A_d = A*dt + eye(n_comb*2);   
B_d = B*dt;
G_d = G*dt;

options = optimoptions('fmincon','Display','iter');
x = fmincon(@fun,K_0,[],[],[],[],[],[],@constr,options)

%Feedback matrix manipulation
K = kron(eye(n_comb),x);
hold on
for i = 1:100

    X_0 = A_d*X_0 + G_d + B_d*K*X_0;
    
    pgon = polyshape(X_0(1,:),X_0(2,:));
    plot(pgon)
    pgon = polyshape(X_0(3,:),X_0(4,:));
    plot(pgon)
    pgon = polyshape(X_0(5,:),X_0(6,:));
    plot(pgon)
    pgon = polyshape(X_0(7,:),X_0(8,:));
    plot(pgon)
end
X_0;
plot(X_0(1,:),X_0(2,:),"*")
plot(X_0(3,:),X_0(4,:),"*")
plot(X_0(5,:),X_0(6,:),"*")
plot(X_0(7,:),X_0(8,:),"*")
plot(X_0(9,:),X_0(10,:),"*")
plot(X_0(11,:),X_0(12,:),"*")
plot(X_0(13,:),X_0(14,:),"*")
plot(X_0(15,:),X_0(16,:),"*")