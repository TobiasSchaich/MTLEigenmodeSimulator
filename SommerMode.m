function [be,ga_0,gc] = SommerMode(a,f,sigma,e2)

if nargin<=3, e2=1;end 
c0 = 299792458;            % m/sec
ep0 = 8.854187817e-12;     % F/m, vacuum permittivity
w = 2*pi*f;
k0 = w/c0;

ec =  1- 1j * sigma./w/ep0;         % sigma must be scalar or size(f)
gc=k0.*sqrt(ec);              % assume good conductor  
J01=besselj(0,gc*a)./besselj(1,gc*a); %calculate exact ratio
J01(isnan(J01))= 1i;            %for large arguments use asymptotic expansion
be = (1-eps)*k0;
ga_0 = sqrt(e2*k0.^2 - be.^2);    
N=1;
while 1          
    gnew = J01.*besselh(1,1,ga_0*a)./besselh(0,1,ga_0*a).*k0.*e2./sqrt(ec);
    if norm(ga_0-gnew)<1e-8, break; end          
    N = N+1;        
    ga_0 = gnew;       
    if N>1e4, break; end         % prevent potential infinite loop
end

ga_0 = gnew;    
gc = sqrt(ga_0.^2 + k0.^2.*(ec-e2));        
be = sqrt(e2*k0.^2 - ga_0.^2);    







