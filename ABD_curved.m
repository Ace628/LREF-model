function [A,B,D,A55]=ABD_curved(z,Neutral_R,b,h,threeDQbar,Layup)

A=0; %A matirx
B=0; %B matrix
D=0; %D matrix
A55=0;

for i=1:size(Layup,1) % The total layers of the laminate

    Qbar_i = squeeze(threeDQbar(i,:,:));
    Qbar_S = [Qbar_i(1,1) Qbar_i(1,2) Qbar_i(1,6)
              Qbar_i(2,1) Qbar_i(2,2) Qbar_i(2,6)
              Qbar_i(6,1) Qbar_i(6,2) Qbar_i(6,6)];

    Ai=Qbar_S.*log((Neutral_R+z(i+1))/(Neutral_R+z(i))); % A matrix
    A=A+Ai;
    
    Bi=Qbar_S.*(z(i+1)-z(i)-Neutral_R*log((Neutral_R+z(i+1))/(Neutral_R+z(i)))); % B matrix
    B=B+Bi;

    Di=Qbar_S.*(0.5*((Neutral_R+z(i+1))^2-(Neutral_R+z(i))^2)-2*Neutral_R*(z(i+1)-z(i))+Neutral_R^2*log((Neutral_R+z(i+1))/(Neutral_R+z(i)))); % D matrix
    D=D+Di;

    A55_i = Qbar_i(5,5).*(-2/h^2*(z(i+1)^2-z(i)^2)+4*Neutral_R/h^2*(z(i+1)-z(i))+(1-4*Neutral_R^2/h^2)*log((Neutral_R+z(i+1))/(Neutral_R+z(i))));
    A55 = A55 + A55_i;

end

A = Neutral_R*b*A;
B = -Neutral_R*b*B;
D = Neutral_R*b*D;
A55 = 5/4*Neutral_R*b*A55;

B(abs(B)<=1e-9)=0;

end
