function [y] = ur_bar_Laminated_REF(x,rot,A,B,D,A55,Neutral_R,K_R,K_R_limit)

% P = 1+(A55*(B(1,1)+Neutral_R*A(1,1))*(D(1,1)+Neutral_R*B(1,1)))/(A(1,1)*D(1,1)-B(1,1)^2)/(B(1,1)-Neutral_R*A55);
P = A55/(B(1,1)-Neutral_R*A55)*(1/A55+(D(1,1)+Neutral_R*B(1,1))/(A(1,1)*D(1,1)-B(1,1)^2)+(Neutral_R^2*A(1,1)+Neutral_R*B(1,1))/(A(1,1)*D(1,1)-B(1,1)^2));
lambda1 = 2-K_R*Neutral_R^2/A55;
lambda2 = 1+K_R*Neutral_R^3*(B(1,1)+Neutral_R*A(1,1))/(A(1,1)*D(1,1)-B(1,1)^2)+...
          K_R*Neutral_R^2*(D(1,1)+Neutral_R*B(1,1))/(A(1,1)*D(1,1)-B(1,1)^2);

X = abs(x - rot);

if X > pi
    X = 2 * pi - X;
end

if K_R >= K_R_limit

alpha = sqrt(-0.5*lambda1+sqrt(0.25*lambda1^2-lambda2));
beta = sqrt(-0.5*lambda1-sqrt(0.25*lambda1^2-lambda2));
Q1 = alpha*(alpha^2+B(1,1)*P+1)/P/(B(1,1)+Neutral_R*A(1,1))-...
     alpha*K_R*Neutral_R^2/P/A55/(B(1,1)+Neutral_R*A(1,1))-...
    (K_R*Neutral_R^3+Neutral_R*A(1,1))/alpha/(B(1,1)+Neutral_R*A(1,1));
Q2 = beta*(beta^2+B(1,1)*P+1)/P/(B(1,1)+Neutral_R*A(1,1))-...
     beta*K_R*Neutral_R^2/P/A55/(B(1,1)+Neutral_R*A(1,1))-...
    (K_R*Neutral_R^3+Neutral_R*A(1,1))/beta/(B(1,1)+Neutral_R*A(1,1));
H1 = (Q1-alpha)/Neutral_R-(alpha^2*Q1+alpha)*(A(1,1)*D(1,1)-B(1,1)^2)/Neutral_R/A55/(D(1,1)+Neutral_R*B(1,1));
H2 = (Q2-beta)/Neutral_R-(beta^2*Q2+beta)*(A(1,1)*D(1,1)-B(1,1)^2)/Neutral_R/A55/(D(1,1)+Neutral_R*B(1,1));

C2 = -Neutral_R*(Q2-Neutral_R*H2)/2/A55/(alpha*Q2-beta*Q1-Neutral_R*alpha*H2+Neutral_R*beta*H1)/sinh(alpha*pi);
C4 = Neutral_R*(Q1-Neutral_R*H1)/2/A55/(alpha*Q2-beta*Q1-Neutral_R*alpha*H2+Neutral_R*beta*H1)/sinh(beta*pi);
C5 = Neutral_R*(B(1,1)+Neutral_R*A(1,1))*(Q1*H2-Q2*H1)/2/pi/A55/(K_R*Neutral_R^2+A(1,1))/(alpha*Q2-beta*Q1-Neutral_R*alpha*H2+Neutral_R*beta*H1);

y = C2*cosh(alpha*X)+C4*cosh(beta*X)+C5;

else

alpha = sqrt(0.5*sqrt(lambda2)-0.25*lambda1);
beta = sqrt(0.5*sqrt(lambda2)+0.25*lambda1);
Q3 = (alpha^2-3*beta^2+B(1,1)*P+1)/P/(B(1,1)+Neutral_R*A(1,1))-K_R*Neutral_R^2/P/A55/(B(1,1)+Neutral_R*A(1,1))-...
     (K_R*Neutral_R^3+Neutral_R*A(1,1))/(alpha^2+beta^2)/(B(1,1)+Neutral_R*A(1,1));
Q4 = (3*alpha^2-beta^2+B(1,1)*P+1)/P/(B(1,1)+Neutral_R*A(1,1))-K_R*Neutral_R^2/P/A55/(B(1,1)+Neutral_R*A(1,1))+...
     (K_R*Neutral_R^3+Neutral_R*A(1,1))/(alpha^2+beta^2)/(B(1,1)+Neutral_R*A(1,1));
H3 = (2*beta^2*Q4-(alpha^2-beta^2)*Q3-1)*(A(1,1)*D(1,1)-B(1,1)^2)/Neutral_R/A55/(D(1,1)+Neutral_R*B(1,1))+(Q3-1)/Neutral_R;
H4 = (2*alpha^2*Q3+(alpha^2-beta^2)*Q4+1)*(A(1,1)*D(1,1)-B(1,1)^2)/Neutral_R/A55/(D(1,1)+Neutral_R*B(1,1))-(Q4-1)/Neutral_R;

C2 = Neutral_R*(alpha*Q3*cosh(alpha*pi)*sin(beta*pi)+beta*Q4*sinh(alpha*pi)*cos(beta*pi))/2/A55/alpha/beta/(Q3-Q4-Neutral_R*H3-Neutral_R*H4)/(cosh(alpha*pi)^2*sin(beta*pi)^2+sinh(alpha*pi)^2*cos(beta*pi)^2)-...
     Neutral_R^2*(alpha*H3*cosh(alpha*pi)*sin(beta*pi)-beta*H4*sinh(alpha*pi)*cos(beta*pi))/2/A55/alpha/beta/(Q3-Q4-Neutral_R*H3-Neutral_R*H4)/(cosh(alpha*pi)^2*sin(beta*pi)^2+sinh(alpha*pi)^2*cos(beta*pi)^2);
C4 = -Neutral_R*(alpha*Q3*sinh(alpha*pi)*cos(beta*pi)-beta*Q4*cosh(alpha*pi)*sin(beta*pi))/2/A55/alpha/beta/(Q3-Q4-Neutral_R*H3-Neutral_R*H4)/(cosh(alpha*pi)^2*sin(beta*pi)^2+sinh(alpha*pi)^2*cos(beta*pi)^2)+...
      Neutral_R^2*(alpha*H3*sinh(alpha*pi)*cos(beta*pi)+beta*H4*cosh(alpha*pi)*sin(beta*pi))/2/A55/alpha/beta/(Q3-Q4-Neutral_R*H3-Neutral_R*H4)/(cosh(alpha*pi)^2*sin(beta*pi)^2+sinh(alpha*pi)^2*cos(beta*pi)^2);
C5 = Neutral_R*(B(1,1)+Neutral_R*A(1,1))*(Q3*H4+Q4*H3)/2/pi/A55/(K_R*Neutral_R^2+A(1,1))/(Q3-Q4-Neutral_R*H3-Neutral_R*H4);

y = C2*cosh(alpha*X)*cos(beta*X)+C4*sinh(alpha*X)*sin(beta*X)+C5;

end

end