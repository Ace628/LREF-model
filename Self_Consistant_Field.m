function [mat] = Self_Consistant_Field(E_f,nu_f,E_m,nu_m,V_f)

% E_f Elastic modulus of the fiber
% E_m Elastic modulus of the matrix
% nu_f Poisson's ratio of the fiber
% nu_m Poisson's ratio of the matrix

G_f = E_f/2/(1+nu_f);
G_m = E_m/2/(1+nu_m);   % shear modulus
K_f = E_f/3/(1-2*nu_f);
K_m = E_m/3/(1-2*nu_m); % bulk modulus

E11 = V_f*E_f+(1-V_f)*E_m;
nu12 = V_f*nu_f+(1-V_f)*nu_m;
nu13 = nu12;
K2 = ((K_f+G_m)*K_m+(K_f-K_m)*G_m*V_f)/((K_f+G_m)-(K_f-K_m)*V_f);
G23 = G_m*(K_m*(G_m+G_f)+2*G_f*G_m+K_m*(G_f-G_m)*V_f)/(K_m*(G_m+G_f)+2*G_f*G_m-(K_m+2*G_m)*(G_f-G_m)*V_f);
E22 = (1/4/K2+1/4/G23+nu12^2/E11)^-1;
E33 = E22;
nu23 = E22/2/G23-1;
G12 = (G_m*(1+V_f)*G_f+(1-V_f)*G_m^2)/((1-V_f)*G_f+(1+V_f)*G_m);
G13 = G12;

mat = [E11 E22 E33 nu12 nu13 nu23 G12 G13 G23];

end