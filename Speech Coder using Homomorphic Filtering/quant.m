function [ quantValue ] = quant(inValue,B, minVal, maxVal)

% 7 bit quantizer on the input value 'inValue' by B bits
% minVal and maxVal denotes the minimum and maximum values of the input
% data :: quantValue
% Reference : http://www.mathworks.com/help/matlab/ref/bsxfun.html
% Reference : https://courses.engr.illinois.edu/ece420/lab3/prelab3.html

% Normalize such that the maxValue of the input comes in [-1 to 1] range
% divide by the smallest power of two such that the resulting absolute value of the largest number
% is less than or equal to one
% This is an easy but fairly reasonable approximation of how numbers 
% outside the range of -1 to 1 are actually handled on the DSP.
N = nextpow2(maxVal-minVal);
inValue = inValue/(2^N);
% Next, quantize to B bits of precision by first multiplying them by 2^{B} rounding to the nearest integer 
Qval = (round((inValue)*(2^(B))))/(2^(B));

% %Set up Quantization levels as a kind of lookup
% lvl = linspace(minVal,maxVal,2^B);
% % Level nearest to the quantized lookup -gives the index from the lookup
% [~,Idx] = min(abs(bsxfun(@minus,inValue,lvl.')));
% quantValue = inValue(Idx);

% Recover orginigal equivalent
quantValue = Qval*(2^N);

end

