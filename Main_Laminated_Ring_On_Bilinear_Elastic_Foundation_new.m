%% If you want to use this code, the citation of our paper is needed
clc
close all
clear all

%% Material properties

%material properties of laminates
%mat(id#,:)=[E11 E22 E33 v12 v13 v23 G12 G13 G23]; % MPa
% mat(1,:)=[20391.63 28.52 28.52 0.39 0.39 0.52 9.64 9.64 9.37];   % Fiber-reinforced rubber
% mat(1,:)=[20384.86 13.06 13.06 0.389 0.389 0.898 3.50 3.50 3.43]; % Fiber-reinforced rubber
% mat(2,:)=[15 15 15 0.48 0.48 0.48 5.068 5.068 5.068]; % RSM rubber
mat(1,:)=[20000 200 200 0.4 0.4 0.6667 100 100 60];

%% Composites
n = 4;         % Ply layers
R = 200;       % The radius of the middle surface of the ring
h = 30;        % The thickness of the whole ring [mm]
b = 60;        % The width of the ring [mm]
Plyt = h/n;    % Ply thickness [mm]
R_in = R-h/2;
Neutral_R = R;

%Layup properties [Material ID, ply thickness, orientation (deg)]
% Layup(1,:)=[1,Plyt,90];
% Layup(2,:)=[1,Plyt,0];
% Layup(3,:)=[1,Plyt,90];
% Layup(4,:)=[1,Plyt,0];
% Layup(5,:)=[1,Plyt,90];

Layup(1,:)=[1,Plyt,0];
Layup(2,:)=[1,Plyt,90];
Layup(3,:)=[1,Plyt,90];
Layup(4,:)=[1,Plyt,0];


%% Neutral axis & ABD matrix
[threeDQbar] = Qbar(Layup,mat);
[z]=mid_z(Layup,h);
[A,B,D,A55]=ABD_curved(z,Neutral_R,b,h,threeDQbar,Layup);

%% Displacements under a point load for the ring on a bilinear elastic foundation
% F = F1(f0); % The magnitude of the concentrated force
K_R = 1.5;       % The radial stiffness of the elastic foundation [N/mm^2]
K_R_C = 0*K_R;  % The compressive stiffness of the elastic foundation [N/mm^2]

% Defined effective stiffnesses and nondimensional parameters
[EA,EV,EI,GA,K_LR,kappa] = eff_stiffness_parameters(A,B,D,A55,Neutral_R,K_R);

N_ang = 2000;          % The number of divided elements
D_ang = 2*pi/N_ang;    % The angle of each element
% theta = [pi/2:D_ang:pi,-pi+D_ang:D_ang:-pi/2];
theta = -pi:D_ang:pi;
theta = theta';
u_r0 = zeros(length(theta),1);
u_theta0 = zeros(length(theta),1);
phi0 = zeros(length(theta),1);
u_r1 = zeros(length(theta),1);
u_theta1 = zeros(length(theta),1);
phi1 = zeros(length(theta),1);
u_r = zeros(length(theta),1);
u_theta = zeros(length(theta),1);
phi = zeros(length(theta),1);
% U_r = zeros(length(theta),1);
% U_theta = zeros(length(theta),1);
% U_phi = zeros(length(theta),1);
rot = 0;    % The rotational degree of coordinate transformation
F_rc = zeros(length(theta),1);
D_F_rc0 = zeros(length(theta),1);
rot_c0 = zeros(length(theta),1);
rot_c1 = zeros(length(theta),1);
r_deform0 = zeros(length(theta),1);
theta_di0 = zeros(length(theta),1);
r_deform1 = zeros(length(theta),1);
theta_di1 = zeros(length(theta),1);

K_R_limit = 4*A55/Neutral_R^2*...
(1+A55*(Neutral_R*B(1,1)+Neutral_R^2*A(1,1))/(A(1,1)*D(1,1)-B(1,1)^2)+A55*(D(1,1)+Neutral_R*B(1,1))/(A(1,1)*D(1,1)-B(1,1)^2));

%% start of the iternation
F1 = 9491.85;
F1 = F1';
U_r1 = zeros(length(F1),1); 
U_r0 = zeros(length(theta),length(F1));
U_theta0 = zeros(length(theta),length(F1));
U_phi0 = zeros(length(theta),length(F1));
K_ring_Nonrecip = zeros(length(F1),1); 
Nor_K_ring_Nonrecip = zeros(length(F1),1); 
Compress_Angle = zeros(length(F1),1);

for f0 = 1:length(F1)
    F = F1(f0);

% Getting the initial compensation load
for x1 = 1:length(theta)

    u_r0(x1) = F*ur_bar_Laminated_REF(theta(x1),rot,A,B,D,A55,Neutral_R,K_R,K_R_limit);
    u_theta0(x1) = F*ut_bar_Laminated_REF(theta(x1),rot,A,B,D,A55,Neutral_R,K_R,K_R_limit);
    phi0(x1) = F*ufi_bar_Laminated_REF(theta(x1),rot,A,B,D,A55,Neutral_R,K_R,K_R_limit);

    theta_di0(x1) = theta(x1) + atan((u_theta0(x1)+(R_in+h-Neutral_R)*phi0(x1))/((R_in+h)+u_r0(x1)));
    r_deform0(x1) = sqrt((R_in+h+u_r0(x1))^2+(u_theta0(x1)+(R_in+h-Neutral_R)*phi0(x1))^2); % Define the initial deformed matrix
   
    if u_r0(x1) >= 0
        D_F_rc0(x1) = 0;
        rot_c0(x1) = 0;
    else
        D_F_rc0(x1) = -(K_R-K_R_C)*u_r0(x1);
        rot_c0(x1) = theta_di0(x1)-sign(theta_di0(x1))*pi;
    end
   
end

U_r = u_r0;
U_theta = u_theta0;
U_phi = phi0;

D_F_rc = D_F_rc0;
rot_c = rot_c0;

theta_di = theta_di0;
r_deform = r_deform0;

% Getting the total deformation
for c0 = 1:5000
    
    D_F_rc1 = D_F_rc;
    F_rc = F_rc + D_F_rc1;
    D_F_rc = zeros(length(theta),1);

    for c1 = 1:length(theta_di)
        for c2 = 1:length(D_F_rc1)
            u_r1(c2) = D_F_rc1(c2)*D_ang*(R_in+h)*ur_bar_Laminated_REF(theta_di(c1),rot_c(c2),A,B,D,A55,Neutral_R,K_R,K_R_limit);
            u_theta1(c2) = D_F_rc1(c2)*D_ang*(R_in+h)*ut_bar_Laminated_REF(theta_di(c1),rot_c(c2),A,B,D,A55,Neutral_R,K_R,K_R_limit);
            phi1(c2) = D_F_rc1(c2)*D_ang*(R_in+h)*ufi_bar_Laminated_REF(theta_di(c1),rot_c(c2),A,B,D,A55,Neutral_R,K_R,K_R_limit);
        end

        u_r(c1) = sum(u_r1);
        u_theta(c1) = sum(u_theta1);
        phi(c1) = sum(phi1);

    end

    U_r = U_r + u_r;
    U_theta = U_theta + u_theta;
    U_phi = U_phi + phi;

    for c3 = 1:length(U_r)
        theta_di0(c3) = theta(c3) + atan((U_theta(c3)+(R_in+h-Neutral_R)*U_phi(c3))/(R_in+h+U_r(c3)));
        r_deform0(c3) = sqrt((R_in+h+U_r(c3))^2+(U_theta(c3)+(R_in+h-Neutral_R)*U_phi(c1))^2); 
    end

    for c4 = 1:length(u_r)
        if u_r(c4) >= 0
           D_F_rc(c4) = 0;
           rot_c1(c4) = 0;
        else
           D_F_rc(c4) = -(K_R-K_R_C)*u_r(c4);
           rot_c1(c4) = theta_di0(c4)-sign(theta_di0(c4))*pi;
        end
    end

% Stop the iteration if the change of the deformation is quite small
    diff = max(abs(u_r));
    if diff < 8e-7
       break
    else
        theta_di1 = theta_di0;
        r_deform1 = r_deform0;

        theta_di = theta_di1;
        r_deform = r_deform1;

        rot_c = rot_c1;
        rot_c1 = zeros(length(theta),1);

        theta_di0 = zeros(length(theta),1);
        r_deform0 = zeros(length(theta),1);

       continue
    end

end

U_r1(f0) = max(abs(U_r));
U_r0(:,f0) = U_r;
K_ring_Nonrecip(f0) = F/U_r1(f0);
Nor_K_ring_Nonrecip(f0) = F/U_r1(f0)/K_LR;
U_theta0(:,f0) = U_theta;
U_phi0(:,f0) = U_phi;

[Min_R, Min_ID] = min(abs(r_deform0 - (R+h/2)));
Compress_Angle(f0) = abs(rad2deg(theta_di0(Min_ID)));

%% Plotting the results
polarplot(theta_di0,r_deform0)

hold on

pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';

end
