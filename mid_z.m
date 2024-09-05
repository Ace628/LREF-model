function [z]=mid_z(Layup,h)
% disp('Calculating the neutral axis')
% Using the midplane of the cross-section is OK to calculate ABD matrix
% The uneven distribution of the stress is considered in B matrix

th=Layup(:,2);  %Layer thicknesses
z=zeros(1,(length(th)+1));

z0 = -h/2;
z(1) = z0;

for j=1:size(Layup,1) % The total layers of the laminate

    z(j+1) = z(j) + th(j);

end

end
