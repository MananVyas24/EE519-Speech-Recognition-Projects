function [ out ] = nonUniformQuant( in )

% Custom non uniform Quantizer
% selects the bit depth of current co-efficient depending upon the dynamic
% range around its neighbours and the average dynamic range

% Since we have the pre-calculated co-efficients we can calculate average
% dynamic range of the whole set

sumDiff = 0;
for i=2:1:length(in)
   diff = abs(in(i) - in(i-1));
   sumDiff = sumDiff + diff;
end
avgDiff = sumDiff/(length(in)-1); % average dynamic range of the whole inp

% Max value of the array of co-effs
maxVal = max(in);
minVal = min(in);

% Based on the Local Dynamic range and the avg dynamic range make a
% decision to fix B
B(1) = 32;
out(1) = quant(in(1),B(1),minVal,maxVal);
for i=2:1:length(in)
   if (abs(in(i)-in(i-1)) > 3*avgDiff)
       B(i) = 64;
       out(i) = quant(in(i),B(i),minVal,maxVal);
   elseif (abs(in(i)-in(i-1)) > 2*avgDiff)
       B(i) = 32;
       out(i) = quant(in(i),B(i),minVal,maxVal);
   elseif (abs(in(i)-in(i-1)) > avgDiff)
       B(i) = 16;
       out(i) = quant(in(i),B(i),minVal,maxVal);
   else
       B(i) = 8;
       out(i) = quant(in(i),B(i),minVal,maxVal);
   end
end

end

