function Sn = generateSine(mag, phase, frq, length)

% % Generated sine wave of the length 'length'
% frq is the normalized frequency,dont pass the frq in Hz !

[m,L]=size(mag);
n=0:(length-1);
s=zeros(1,length);

for k=1:L
    temp=mag(k)*cos((frq(k)*n)+phase(k));
    s=s+temp;
end
Sn=s;
end % of function