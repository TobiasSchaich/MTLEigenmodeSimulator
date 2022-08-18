function [T_I,be,Zc,V_modal,I_modal] = CalculateEigenmodes(R,C,L,G,f)
%CalculateEigenmodes Determines Eigenmodes in current and Voltage and
%calculates propagation constant (e^-1i*be*z) and characteristic impedance
%   R,C,L,G - TL parameters, 
%   f  - Frequency in Hz
%   T  - matrix w normalised modal currents as column vectors/Transformation-matrix
%   be - propagation constants matrix form
%   Zc - characteristic impedance matrix relating modal voltages to modal currents
%   V_modal - modal voltages normalised to 1 in terminal basis
%   I_modal - modal currents in terminal basis
        w=2*pi*f;
        Z=R+1i*w*L;
        Y=G+1i*w*C;
        [T_I,ga_square]=eig(Y*Z); 
        %[T_V,~]=eig(Z*Y); 
        ga=sqrt(ga_square); % d/dz I=ga^2*I --> exp(-ga*z) for modal currents; neglect backwards propagating waves
        [~, idx]=sort(diag(ga)); %Sort by losses / lowest loss first
        ga=ga(idx,idx); 
        T_I=T_I(:,idx);
        T_V=inv(T_I).';
        Zc=(T_I.')*Z*T_I*inv(ga); % T_I.'=inv(T_V), Relation: V_modal=Z_c*I_modal
        be=-1i*ga; %exp(ga*z)=exp(-i*be*z) -->be=-1i*ga
        %bets=diag(be); k0=2*pi*1e9/299792458;
        %dat=[real(bets(1))/k0, real(bets(2))/k0, imag(bets(1))*20*log10(exp(1)), imag(bets(2))*20*log10(exp(1))]
        %dlmwrite('Data.csv',dat,'delimiter',',','precision',10, '-append')
        V_modal=T_V;
        I_modal=T_I; 
end

