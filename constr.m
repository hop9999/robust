function [c,ceq] = constr(x)
%Modeling params
dt = 1e-2;
N = 100;

%trajectory constraints
x_lim = 0.15;
dx_lim = 0.25;
u_lim = 6;
x_f = 0.05;
dx_f = 0.05;

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

%Feedback matrix manipulation
K = kron(eye(n_comb),x);
 
%init coordinates
N_init_x = 4;
x_0 = [0.1 -0.1 -0.1 0.1 
       -0.1 -0.1 0.1 0.1];
   
X_0 = kron(ones(n_comb,1),x_0);

%Dynamic discretization
A_d = A*dt + eye(n_comb*2);   
B_d = B*dt;
G_d = G*dt;

%constraints preparing
c = zeros(((256 + n_comb*N_init_x*2)*N + 256),1);
ceq = [];
A_lim_x = [eye(64) zeros(64,64)
           -eye(64) zeros(64,64)
           zeros(64,64) eye(64)
           zeros(64,64) -eye(64)];
A_lim_u = [eye(16*4);-eye(16*4)];

for i = 1:N
    
    U = K*X_0;
    X_0 = A_d*X_0 + G_d + B_d*U;
    
    X = reshape(X_0,[2,64]);
    X = reshape(X',[1,128])';
    U = reshape(U,[1,16*4])';
    
    c((256 + n_comb*N_init_x*2)*(i-1) + 1:(256 + n_comb*N_init_x*2)*i,1) = [A_lim_x*X - [x_lim*ones(128,1);dx_lim*ones(128,1)];
                                                                            A_lim_u*U - u_lim*ones(n_comb*N_init_x*2,1)];

end
c(end - 255:end) = [A_lim_x*X - [x_f*ones(128,1);dx_f*ones(128,1)]];
