function [R,C,L,G] = CalculateRCLG_coated_NoGND(f, pos_vec, a, b,  ed_diel,  sigma_cond,LossTanDiel)
%UNTITLED2 Calculates new RCLG parameters
%   Detailed explanation goes here
mu0=4*pi*1e-7;
eps0=8.854187*1e-12;

N_cond=length(pos_vec(1,:)); %number of conductors in system excluding ground
k0=2*pi*f*sqrt(mu0*eps0);
[be] = GoubauMode(a,b,ed_diel,f);
ga=sqrt(k0^2-be^2); 
gd=sqrt(k0^2*ed_diel-be^2); 
%Resistance
skin_depth=sqrt(1./(sigma_cond*pi*f*mu0)); 
R=eye(N_cond)*1./(2*pi*a*sigma_cond*skin_depth);
% Capacitance/Inductance
Z0=@(r) besselj(0,gd*r).*bessely(0,gd*a)-bessely(0, gd*r).*besselj(0,gd*a);
Z1=@(r) besselj(1,gd*r).*bessely(0,gd*a)-bessely(1, gd*r).*besselj(0,gd*a);

P_self=eye(N_cond).*(-1/(2*pi*eps0*ed_diel*(1-1i*LossTanDiel)*gd*a).*Z0(b)./Z1(a)-1/(2*pi*eps0*ga*a)*...
    Z1(b)/Z1(a)*(-besselh(0,ga*b))/besselh(1,ga*b));
L_self=eye(N_cond).*(-mu0/(2*pi*gd*a).*Z0(b)./Z1(a)-mu0/(2*pi*ga*a)*...
    Z1(b)/Z1(a)*(-besselh(0,ga*b))/besselh(1,ga*b));

P_mutual=zeros(N_cond);
L_mutual=zeros(N_cond); 
for ind1=1:N_cond
    for ind2=(ind1+1):N_cond
        P_mutual(ind1,ind2)=-1/(2*pi*a*ga*eps0)*Z1(b)/Z1(a).*...
            (-besselh(0,ga*norm(pos_vec(:,ind1)-pos_vec(:,ind2))))/besselh(1,ga*b);
        L_mutual(ind1,ind2)=-mu0/(2*pi*ga*a)*Z1(b)/Z1(a)*...
            (-besselh(0,ga*norm(pos_vec(:,ind1)-pos_vec(:,ind2))))/besselh(1,ga*b);
    end
end


P=P_self+P_mutual+P_mutual.';
C_complex=inv(P);
G=-2*pi*f*imag(C_complex); 
C=real(C_complex);
L=L_self+L_mutual+L_mutual.';
end
