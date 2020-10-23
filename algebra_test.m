clear all;
clc;
close all;

dt = 1e-2;

k_i = 0.03;
mgl = 0.183;
I = 0.0040;
k_q = 0.001;

x_d = pi/2;

A_c = [0                1
       -mgl*cos(x_d)/I  -5*k_q/I];

A_d1 = [0  0
        0  -5*k_q/I];
   
A_d2 = [0                    0
        -0.1*mgl*cos(x_d)/I  0];

B_c = [0 
       k_i/I];

B_d = [0 
       0.1*k_i/I];
 
G_c = [0
       0];
     
G_d = [0
       0.1*mgl*sin(x_d)/I];
     
n_comb = 16;

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

K = [-1 -1];  

K = kron(eye(n_comb),K)

x_0 = [1 -1 -1 1 
       -0.1 -0.1 0.1 0.1];
   
X_0 = kron(ones(n_comb,1),x_0)
   

hold on
pgon = polyshape(X_0(1,:),X_0(2,:));
plot(pgon)

A_d = A*dt + eye(n_comb*2);
          
B_d = B*dt;
G_d = G*dt;

for i = 1:1
    U = K*X_0
    X_0 = A_d*X_0 + G_d + B_d*U;
    U = reshape(U,[1,16*4])';
    A = [eye(16*4);-eye(16*4)];
    A*U - 1*ones(16*8,1);
    X_k = [X_0(1,:),X_0(3,:),X_0(5,:),X_0(7,:),...
           X_0(9,:),X_0(11,:),X_0(13,:),X_0(15,:),...
           X_0(17,:),X_0(19,:),X_0(21,:),X_0(23,:),...
           X_0(25,:),X_0(27,:),X_0(29,:),X_0(31,:);...
           X_0(2,:),X_0(4,:),X_0(6,:),X_0(8,:),...
           X_0(10,:),X_0(12,:),X_0(14,:),X_0(16,:),...
           X_0(18,:),X_0(20,:),X_0(22,:),X_0(24,:),...
           X_0(26,:),X_0(28,:),X_0(30,:),X_0(32,:)];
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