function [be_vector,ga_0,ga_d] = GoubauMode(a,b,ed,f)

c0 = 299792458;     
w = 2*pi*f;
k0 = w/c0;           
kd = k0*sqrt(ed);

gd = @(be) sqrt(kd.^2 - be.^2);
h = @(be) sqrt(be.^2 - k0.^2); 


Z0 = @(be) besselj(0, gd(be).*b).*bessely(0, gd(be).*a)-...
    bessely(0, gd(be).*b).*besselj(0, gd(be).*a);
Z1 = @(be) besselj(1, gd(be).*b).*bessely(0, gd(be).*a)-...
    bessely(1, gd(be).*b).*besselj(0, gd(be).*a);

fun_vector = @(be) (ed./gd(be).*Z1(be)./Z0(be) + ...
    1./h(be).*besselk(1, h(be).*b)./besselk(0, h(be).*b)).'; %function
n=(1+eps);
%be_interval=(1-eps).*kd.'
if length(b)==1
    be_vector=zeros(1,length(f));
    for idx=1:length(f)
        fun_scalar=@(be) subsref(fun_vector(be), struct('type', '()', 'subs', {{idx}}));
        be_interval=[n*k0(idx),(1-eps).*kd(idx)]; 
        be=fzero(fun_scalar, be_interval);
        n=be./k0(idx);
        be_vector(idx)=be; 
    end
    ga_0 = sqrt(k0.^2 - be_vector.^2);
    ga_d = gd(be_vector);
elseif length(f)==1
    be_vector=zeros(1,length(b));
    for idx=1:length(b)
        be_interval=[n*k0,(1-eps).*kd]; 
        fun_scalar=@(be) subsref(fun_vector(be), struct('type', '()', 'subs', {{idx}}));
        be=fzero(fun_scalar, be_interval);
        n=be./k0;
        be_vector(idx)=be; 
    end
    ga_0 = sqrt(k0.^2 - be_vector.^2);
    ga_d = gd(be_vector);
end
















