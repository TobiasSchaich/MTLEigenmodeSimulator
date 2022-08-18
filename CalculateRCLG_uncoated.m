function [R,C,L,G] = CalculateRCLG_uncoated(f, pos_vec, rad_cond,  eps_eff_rel, sigma_diel, sigma_cond)
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
%L_self=real(diag(-mu0/(2*pi*ga*rad_cond)*(besselh(0,2*ga*pos_vec(2,:))-besselh(0,ga*rad_cond))/besselh(1,ga*rad_cond))); % self inductance
L_self=-mu0/2/pi/rad_cond*real(1i/ga^2*(besselh(0,2*ga*pos_vec(2,:))-besselh(0,ga*rad_cond)))/real(1i/ga*besselh(1,ga*rad_cond)); 
pos_image_vec=CreateImageVec(pos_vec);
for ind1=1:N_cond
    for ind2=(ind1+1):N_cond
        L_mutual(ind1,ind2)=real(-mu0/(2*pi*ga*rad_cond)*...
            (besselh(0,ga*norm(pos_image_vec(:,ind1)-pos_vec(:,ind2)))-besselh(0,ga*norm(pos_vec(:,ind1)-pos_vec(:,ind2))))/besselh(1,ga*rad_cond));
    end
end

L=L_self+L_mutual+L_mutual.';
P=1/(eps0*eps_eff_rel*mu0)*L;
C=inv(P);
G=mu0*sigma_diel*inv(L); 
delta=sqrt(1/(pi*f*mu0*sigma_cond));
R=eye(N_cond)*1/(2*pi*rad_cond*delta*sigma_cond);
%L=L+eye(N_cond)*1/(4*pi*rad_cond)*sqrt(mu0/(pi*sigma_cond*f));
end



function pos_image_vec=CreateImageVec(pos_vec)
    pos_image_vec=pos_vec;
    pos_image_vec(2,:)=-pos_image_vec(2,:);
end