function [threeDQbar] = Qbar(Layup,mat)

%mat(id#,:)=[E11 E22 E33 v12 v13 v23 G12 G13 G23]; % MPa
%Layup properties [Material ID, ply thickness(mm), orientation(deg)]

Q = zeros(6,6); % Q matrix
R=eye(6);
R(4,4)=2;
R(5,5)=2;
R(6,6)=2; % R matrix
e=10^(-9); % error threshold

numArrays = size(Layup,1);
threeDQbar = zeros(numArrays,6,6);

for i=1:size(Layup,1) % The total layers of the laminate
    M_id=Layup(i,1); % Material ID

    % Material properties
    E11 = mat(M_id,1);
    E22 = mat(M_id,2);
    E33 = mat(M_id,3);
    v12 = mat(M_id,4);
    v13 = mat(M_id,5);
    v23 = mat(M_id,6);
    G12 = mat(M_id,7);
    G13 = mat(M_id,8);
    G23 = mat(M_id,9);
    v21 = E22*v12/E11;
    v31 = E33*v13/E11;
    v32 = E33*v23/E22;

    % Q matrix
    dum = 1-v12*v21-v23*v32-v31*v13-2*v21*v32*v13; %1-v12*v21=1-v12^2*E22/E11 v12*E22/E11=v21
    Q(1,1) = E11*(1-v23*v32)/dum;
    Q(2,2) = E22*(1-v31*v13)/dum;
    Q(3,3) = E33*(1-v12*v21)/dum;
    Q(4,4) = G23;
    Q(5,5) = G13;
    Q(6,6) = G12;
    Q(1,2) = (v21+v31*v23)*E11/dum;
    Q(2,1) = Q(1,2);
    Q(1,3) = (v31+v21*v32)*E11/dum;
    Q(3,1) = Q(1,3);
    Q(2,3) = (v32+v12*v31)*E22/dum;
    Q(3,2) = Q(2,3);

    % T matrix
    angle=pi*Layup(i,3)/180; %Unit transfer: deg->rad
    m=cos(angle);
    n=sin(angle);
    T = [m^2 n^2 0 0 0 2*m*n
         n^2 m^2 0 0 0 -2*m*n
         0 0 1 0 0 0
         0 0 0 m -n 0
         0 0 0 n m 0
         -m*n m*n 0 0 0 m^2-n^2]; % Transformation matrix

    Qb = T\Q*R*T/R;

    threeDQbar(i,:,:) = Qb;
end

end