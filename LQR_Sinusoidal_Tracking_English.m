%Initialization 
clear all; close all; clc;
%System configuration 
m_sys = 1; b_sys = 0.5; k_sys = 1;
%Establish system matrix 
A = [0 1; -k_sys/m_sys -b_sys/m_sys];
n = size(A,1);
B = [0;1/m_sys];
p = size(B,2);
C = [1 0];
D = [0];

%Time step 
Ts = 0.01;
%Discretize the system 
sys_d = c2d(ss(A,B,C,D), Ts);
%Get discrete system A and B matrix 
A = sys_d.a; B = sys_d.b;

%Weight matrices design 
%State 
Q = [10000000 0;0 0.5];
S = [10000000 0;0 0.5];
%Input 
R = 0.08;

%Define number of steps 
n_steps = 100;
%Define control goal 
w = 2*pi;
A_goal = [0 w;-w 0]; B_goal = [0;1];%State space form of sinusoidal signal given unit impulse input 
C_goal = [1 0]; D_goal = [0];
%System discretization 
sys_goal_d = c2d(ss(A_goal,B_goal,C_goal,D_goal),Ts);
%Get discrete system A and B matrix 
A_goal_d = sys_goal_d.a; B_goal_d = sys_goal_d.b;
xm_goal = zeros(2,n_steps); um_goal = (1/Ts)*[1 zeros(1,n_steps-1)];%Establish a unit impulse function(The reason of multiplying the impulse signal 
% with a (1/Ts) is to correct the output amplitude based on observed output
% amplitude and desired trajectory amplitude)

for j = 2:n_steps
    xm_goal(:,j) = A_goal_d*xm_goal(:,j-1) + B_goal_d*um_goal(j-1);
end

%Get desired input at every time step and store all the values in variable ud 
ud = zeros(p,n_steps);
for k = 1:n_steps
    ud(k) = mldivide(B,(A_goal_d-A)*xm_goal(:,k));
end

x_history = zeros(n,n_steps);
u_history = zeros(p,n_steps);

%Prepare augmented matrices 
Ca = [eye(n) -eye(n)];Aa=[A A_goal_d-A;zeros(n) A_goal_d];
Ba = [B;zeros(n,p)];
Qa = transpose(Ca)*Q*Ca;
Sa = transpose(Ca)*S*Ca;
%Feedback gain calculation 
for k = 1:n_steps
 F = inv(R+Ba'*Sa*Ba)*Ba'*Sa*Aa;
 Sa = (Aa-Ba*F)'*Sa*(Aa-Ba*F)+F'*R*F+Qa;
 if k == 1
     F_N = F;
 else
     F_N = [F;F_N];
 end
end

x = [0;0];u = 0;
%Simulation 
for k = 1:n_steps
    xa = [x;xm_goal(:,k)];
    du = -F_N(k,:)*xa;
    u = du + ud(k);
    x = A*x + B*u;
    x_history(:,k) = x;
    u_history(:,k) = u;
end
%Plot results, only displacement is drawn in this case 
figure(1);clf();plot(1:n_steps,x_history(1,:),'-r','LineWidth',1.5);hold on;
plot(1:n_steps,xm_goal(1,:),'--b','LineWidth',1.5);
xlabel('Iteration Steps');ylabel('Displacement(m)');legend('System Response','Control Goal');
title('System Displacement Response');