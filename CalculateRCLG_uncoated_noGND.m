function [R,C,L,G] = CalculateRCLG_uncoated_noGND(f,pos_vec, rad_cond, eps_eff_rel, sigma_diel, sigma_cond)
%CalculateRCLG_Old calculates RCLG using electrostatic approach
%   pos_vec - positions of conductors relative to ground
%   rad_cond - radius of conductors
%   eps_eff_rel - effective dielectric constnat
%   sigma_diel - conductivity of dielectric for incorporation of losses
%   sigma_cond - conductivity of conductor for incorporation of losses
%   f - frequency
mu0=4*pi*1e-7;% H/m
eps0=8.85418781*1e-12;%F/m
N_cond=length(pos_vec(1,:)); %number of conductors in system excluding ground
L=zeros(N_cond); L_mutual=L; 
[~,ga,~]=SommerMode(rad_cond,f,sigma_cond,eps_eff_rel);
L_self=eye(N_cond).*(mu0/(2*pi*rad_cond))*(besselh(0,ga*rad_cond)/(ga*besselh(1,ga*rad_cond))); % self inductance
for ind1=1:N_cond
    for ind2=(ind1+1):N_cond
        L_mutual(ind1,ind2)=(mu0/(2*pi*rad_cond))*...
        (1/ga*besselh(0,ga*norm(pos_vec(:,ind1)-pos_vec(:,ind2))))/(besselh(1,ga*rad_cond));
    end
end

L=L_self+L_mutual+L_mutual.';

P=1/(eps0*eps_eff_rel*mu0)*L;
C=inv(P);
G=mu0*sigma_diel.*inv(L); 
delta=sqrt(1/(pi*f*mu0*sigma_cond));
R=eye(N_cond)*(1+1i)./(2*pi*rad_cond*delta*sigma_cond);
end
