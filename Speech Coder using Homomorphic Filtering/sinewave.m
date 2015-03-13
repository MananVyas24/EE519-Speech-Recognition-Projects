function result=sinewave(A,f,phase,length)
%
%   result = synthesized speech using sinewaves
%   A      = measured magnitudes
%   f      = normalized freq.
%   phase  = measured phase
%   length = duration of the synthesized speech

[m,L]=size(A);
n=0:(length-1);
s=zeros(1,length);

for k=1:L
    temp=A(k)*cos((f(k)*n)+phase(k));
    s=s+temp;
end

result=s;