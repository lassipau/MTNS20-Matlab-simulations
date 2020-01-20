function [G1,G2] = ConstrIMReal(freqs,dimY)
% G1 = ConstrIMReal(freqs,dimY)
%
% Construct a real block-diagonal internal model
% freqs = Frequencies to be included in the controller, only real nonnegative
% frequencies, if zero frequency is included, it's the first element in the
% vector
% dimY = Number of copies of each frequency to be included in G1
% G1 = System matrix of the internal model
%
% Copyright (C) 2020 by Lassi Paunonen (lassi.paunonen@tuni.fi)
% Licensed under GNU GPLv3 (see LICENSE.txt).

q = length(freqs);

if freqs(1)==0, dimZ = dimY*(2*q-1); else, dimZ = dimY*2*q; end  
  
G1 = zeros(dimZ);

if freqs(1)==0
  zoffset = dimY; 
  nzfreqs = freqs(2:end);
  
  G2 = [eye(dimY);repmat([eye(dimY);zeros(dimY)],q-1,1)];

else
  zoffset = 0;
  nzfreqs = freqs;
  
  G2 = repmat([eye(dimY);zeros(dimY)],q,1);

end

for ind = 1:length(nzfreqs)
  indran = zoffset+(ind-1)*2*dimY+(1:(2*dimY));

  G1(indran,indran) = nzfreqs(ind)*[zeros(dimY) eye(dimY);-eye(dimY) zeros(dimY)];
  
end

