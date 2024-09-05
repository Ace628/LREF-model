clc
close all
clear all

%% Material properties
% Material properties of fibers and matrix
E_f = 40000;   % Elastic modulus of fibers [MPa]
nu_f = 0.22;   % Poisson's ratio of fibers
E_m = 5.0;    % Elastic modulus of matrix [MPa] eta=1 E_m=4.834
nu_m = 0.48;   % Poisson's ratio of matrix
V_f = 0.5;     % The volume fraction of fibers

%material properties of laminates
%mat(id#,:)=[E11 E22 E33 v12 v13 v23 G12 G13 G23]; % MPa
% mat(1,:)=[20391.63 28.52 28.52 0.39 0.39 0.52 9.64 9.64 9.37]; % Fiber-reinforced rubber
% mat(1,:)=[20384.86 13.06 13.06 0.389 0.389 0.898 3.50 3.50 3.43];
% mat(1,:)=Self_Consistant_Field(E_f,nu_f,E_m,nu_m,V_f);
% mat(1,:)=[6740 1142.61 1142.61 0.45 0.45 0.48 386.23 386.23 385.80];
mat(1,:)=[20000 20 20 0.4 0.4 0.6667 10 10 6];
% mat(1,:)=[20000 20000 20000 0.4 0.4 0.4 7142.857 7142.857 7142.857];
mat(2,:)=[15 15 15 0.48 0.48 0.48 5.068 5.068 5.068]; % RSM rubber

%% Composites
n = 4;    % Ply layers
R = 200;    % The radius of the middle surface of the ring
h = 20;      % The thickness of the whole ring [mm]
R_in = R-h/2;
b = 60;       % The width of the ring [mm]
Plyt = h/n;    % Ply thickness [mm]
Neutral_R = R;

%Layup properties [Material ID, ply thickness, orientation (deg)]
Layup(1,:)=[1,Plyt,0];
Layup(2,:)=[1,Plyt,90];
Layup(3,:)=[1,Plyt,90];
Layup(4,:)=[1,Plyt,0];
% Layup(5,:)=[1,Plyt,90];

%% ABD matrix
[threeDQbar] = Qbar(Layup,mat);
[z]=mid_z(Layup,h);
[A,B,D,A55]=ABD_curved(z,Neutral_R,b,h,threeDQbar,Layup);

%% Displacements under opposite concentrated forces
% K_R1 = 0.1:0.1:15;
K_R1 = 20;
K_ring_k = zeros(length(K_R1),1);
Nor_K_ring_k = zeros(length(K_R1),1);
kappa1 = zeros(length(K_R1),1);
for k0 = 1:length(K_R1)
K_R = K_R1(k0);    % The radial stiffness of the elastic foundation [N/mm^2]

% Defined effective stiffnesses and nondimensional parameters
[EA,EV,EI,GA,K_LR,kappa] = eff_stiffness_parameters(A,B,D,A55,Neutral_R,K_R);

theta = -pi:pi/1000:pi;
theta = theta';
u_r = zeros(length(theta),1);
u_theta = zeros(length(theta),1);
phi = zeros(length(theta),1);
theta_di = zeros(length(theta),1);
r_deform = zeros(length(theta),1);
rot = 0;    % The rotational degree of coordinate transformation

K_R_limit = 4*A55/Neutral_R^2*...
(1+A55*(Neutral_R*B(1,1)+Neutral_R^2*A(1,1))/(A(1,1)*D(1,1)-B(1,1)^2)+A55*(D(1,1)+Neutral_R*B(1,1))/(A(1,1)*D(1,1)-B(1,1)^2));

% P1 = 1+(A55*(B(1,1)+Neutral_R*A(1,1))*(D(1,1)+Neutral_R*B(1,1)))/(A(1,1)*D(1,1)-B(1,1)^2)/(B(1,1)-Neutral_R*A55);
% P2 = B(1,1)*A55/(B(1,1)-Neutral_R*A55)*(1/A55+(D(1,1)+2*Neutral_R*B(1,1)+Neutral_R^2*A(1,1))/(A(1,1)*D(1,1)-B(1,1)^2));
% dum0 = A(1,1)*D(1,1)-B(1,1)^2;

F1 = 10000;
F1 = F1';

U_r0 = zeros(length(theta),length(F1));
Deflection = zeros(length(F1),1);
K_ring = zeros(length(F1),1);
Nor_K_ring = zeros(length(F1),1);
U_th0 = zeros(length(theta),length(F1));
U_phi0 = zeros(length(theta),length(F1));

for f0 = 1:length(F1)
F = F1(f0); % The magnitude of the concentrated force

for x = 1:length(theta)
      
    u_r(x) = F*ur_bar_Laminated_REF(theta(x),rot,A,B,D,A55,Neutral_R,K_R,K_R_limit);
    u_theta(x) = F*ut_bar_Laminated_REF(theta(x),rot,A,B,D,A55,Neutral_R,K_R,K_R_limit);
    phi(x) = F*ufi_bar_Laminated_REF(theta(x),rot,A,B,D,A55,Neutral_R,K_R,K_R_limit);

    theta_di(x) = theta(x) + atan((u_theta(x)+(R_in+h-Neutral_R)*phi(x))/((R_in+h)+u_r(x)));
    r_deform(x) = sqrt((R_in+h+u_r(x))^2+(u_theta(x)+(R_in+h-Neutral_R)*phi(x))^2);

end

U_r0(:,f0) = u_r;
Deflection(f0) = max(abs(u_r));
U_th0(:,f0) = u_theta;
U_phi0(:,f0) = phi;
K_ring(f0) = F/Deflection(f0);
Nor_K_ring(f0) = K_ring(f0)/K_LR;

figure(1)
plot(theta,u_r)
figure(2)
% plot(theta,u_theta)
% figure(3)
% plot(theta,phi)
% figure(4)
polarplot(theta_di,r_deform)
hold on
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';

end

% kappa1(k0) = kappa;
% K_ring_k(k0) = K_ring;
% Nor_K_ring_k(k0) = Nor_K_ring;

end
% figure(1)
% plot(kappa1,K_ring_k)
% figure(2)
% plot(kappa1,Nor_K_ring_k)

% figure(2)
% plot(Deflection,F1);