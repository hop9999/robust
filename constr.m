function [c,ceq] = constr(x)
dt = 1e-2;
x_lim = 0.1;
dx_lim = 1.2;
u_lim = 6;
x_f = 0.05;
dx_f = 0.05;

k_i = 0.031;
mgl = 0.1852;
I = 0.0046;
k_q = 0.0012;
x_d = pi/2;

A = [0  1 0  0
    -mgl*cos(x_d)/I -k_q/I 0 0
    0  0 0  1
    0  0 -mgl*cos(x_d)/I -k_q/I];

B = [0  0 
     k_i/I 0
    0  0 
     0 k_i/I];

G = 0.12*[0
     mgl*sin(x_d)/I
     0
     -mgl*sin(x_d)/I];
 


X_0 = 0.1*[2 -2 -2 2 
           -1 -1 1 1
           2 -2 -2 2 
           -1 -1 1 1];

K = [x(1) x(2) 0 0
     0 0 x(1) x(2)];


A_d = A*dt + [eye(2),zeros(2,2);zeros(2,2),eye(2)];
B_d = B*dt;
G_d = G*dt;

c = [];
ceq = [];
N_v = 20;
for i = 1:200

    X_0 = A_d*X_0 + G_d + B_d*K*X_0;
    %X_0 = reshape(X_0,[2,length(X_0)*2])';
    %k = convhull(X_0);
    %X_0 = [X_0(k,1),X_0(k,2)]';
    if length(X_0) > N_v
         X_n1 = interp1( 1:1:length(X_0),X_0(1,:),linspace(1,length(X_0),N_v));
         X_n2 = interp1( 1:1:length(X_0),X_0(2,:),linspace(1,length(X_0),N_v));
         X_n3 = interp1( 1:1:length(X_0),X_0(3,:),linspace(1,length(X_0),N_v));
         X_n4 = interp1( 1:1:length(X_0),X_0(4,:),linspace(1,length(X_0),N_v));
         X_0 = [X_n1;X_n2;X_n3;X_n4];
    end
    X_k = [X_0(1,:),X_0(3,:);X_0(2,:),X_0(4,:)];
    c = [c, max(X_k(1,:))- x_lim, -min(X_k(1,:))- x_lim, max(X_k(2,:))- dx_lim, -min(X_k(2,:))- dx_lim, max(x*X_k) - u_lim, -min(x*X_k) - u_lim];
    
end
c = [max(X_k(1,:))- x_f, -min(X_k(1,:))- x_f, max(X_k(2,:))- dx_f, -min(X_k(2,:))- dx_f];
