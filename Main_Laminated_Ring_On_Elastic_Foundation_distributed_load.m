clc
close all
clear all

%% Material properties
% Material properties of fibers and matrix
E_f = 40000;   % Elastic modulus of fibers [MPa]
nu_f = 0.22;   % Poisson's ratio of fibers
E_m = 5.0;     % Elastic modulus of matrix [MPa] eta=1 E_m=4.834
nu_m = 0.48;   % Poisson's ratio of matrix
V_f = 0.5;     % The volume fraction of fibers

%material properties of laminates
%mat(id#,:)=[E11 E22 E33 v12 v13 v23 G12 G13 G23]; % MPa
% mat(1,:)=[20391.63 28.52 28.52 0.39 0.39 0.52 9.64 9.64 9.37]; % Fiber-reinforced rubber
% mat(1,:)=[20384.86 13.06 13.06 0.389 0.389 0.898 3.50 3.50 3.43];
% mat(1,:)=Self_Consistant_Field(E_f,nu_f,E_m,nu_m,V_f);
% mat(1,:)=[6740 1142.61 1142.61 0.45 0.45 0.48 386.23 386.23 385.80];
% mat(1,:)=[31000 8619.71 8619.71 0.35 0.35 0.58 2971.58 2971.58 2721.41];
% mat(2,:)=[15 15 15 0.48 0.48 0.48 5.068 5.068 5.068]; % RSM rubber
mat(1,:)=[20000 200 200 0.4 0.4 0.6667 100 100 60];

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

%% ABD matrix
[threeDQbar] = Qbar(Layup,mat);
[z]=mid_z(Layup,h);
[A,B,D,A55]=ABD_curved(z,Neutral_R,b,h,threeDQbar,Layup);

K_R = 0.1;    % The radial stiffness of the elastic foundation [N/mm^2]
F = 60; 

% Defined effective stiffnesses and nondimensional parameters
[~,~,~,~,K_LR,kappa] = eff_stiffness_parameters(A,B,D,A55,Neutral_R,K_R);
theta = -pi:pi/1000:pi;
theta = theta';
u_r1 = zeros(length(theta),1);
u_theta1 = zeros(length(theta),1);
phi1 = zeros(length(theta),1);
u_r = zeros(length(theta),1);
u_theta = zeros(length(theta),1);
phi = zeros(length(theta),1);
theta_di = zeros(length(theta),1);
r_deform = zeros(length(theta),1);

rot_F = deg2rad(30);
N_ang = 1000;           % The number of divided elements in the loading region
D_ang = 2*rot_F/N_ang;
rot0 = -rot_F:D_ang:rot_F;
rot0 = rot0';  % The rotational degree of coordinate transformation
% D_F0 = F*cos(rot0*pi/2/max(rot0));         % Case 1 Cosine function
D_F0 = F*abs(sin(rot0*pi/2/max(rot0)));  % Case 2 Sine function
% D_F0 = F*ones(length(rot0),1);

K_R_limit = 4*A55/Neutral_R^2*...
(1+A55*(Neutral_R*B(1,1)+Neutral_R^2*A(1,1))/(A(1,1)*D(1,1)-B(1,1)^2)+A55*(D(1,1)+Neutral_R*B(1,1))/(A(1,1)*D(1,1)-B(1,1)^2));

%% Getting the deformation under the distributed load
for c1 = 1:length(theta)
    for c2 = 1:length(D_F0)
    u_r1(c2) = D_F0(c2)*D_ang*(R_in+h)*ur_bar_Laminated_REF(theta(c1),rot0(c2),A,B,D,A55,Neutral_R,K_R,K_R_limit);
    u_theta1(c2) = D_F0(c2)*D_ang*(R_in+h)*ut_bar_Laminated_REF(theta(c1),rot0(c2),A,B,D,A55,Neutral_R,K_R,K_R_limit);
    phi1(c2) = D_F0(c2)*D_ang*(R_in+h)*ufi_bar_Laminated_REF(theta(c1),rot0(c2),A,B,D,A55,Neutral_R,K_R,K_R_limit);
    end
    u_r(c1) = sum(u_r1);
    u_theta(c1) = sum(u_theta1);
    phi(c1) = sum(phi1);
end

for c3 = 1:length(u_r)
    theta_di(c3) = theta(c3) + atan((u_theta(c3)+(R_in+h-Neutral_R)*phi(c3))/((R_in+h)+u_r(c3)));
    r_deform(c3) = sqrt((R_in+h+u_r(c3))^2+(u_theta(c3)+(R_in+h-Neutral_R)*phi(c3))^2);
end

% figure(1)
% plot(theta,u_r)
% figure(2)
% plot(theta,u_theta)
% figure(3)
% plot(theta,phi)
% figure(4)
polarplot(theta_di,r_deform)
% hold on
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
