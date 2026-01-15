function [EA,EV,EI,GA,K_LR,kappa] = eff_stiffness_parameters(A,B,D,A55,Neutral_R,K_R)

EA = A(1,1)-B(1,1)^2/D(1,1);
EV = A(1,1)*D(1,1)/B(1,1)-B(1,1);
EI = D(1,1)-B(1,1)^2/A(1,1);
GA = A55;

K_LR = 8*pi/(pi^2*Neutral_R/GA+pi^2*(Neutral_R/EA+Neutral_R^2/EV)+(pi^2-8-8*EI/Neutral_R/EV)*(Neutral_R^3/EI+Neutral_R^2/EV));
kappa = Neutral_R*K_R/K_LR;


end
