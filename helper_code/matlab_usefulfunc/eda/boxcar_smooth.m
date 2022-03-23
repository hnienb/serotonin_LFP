function newvec = boxcar_smooth(oldvec, carlen)
% apply boxcar convolution to smooth the given vector
% this function automatically removes edge noise

car = ones(1, carlen)./carlen;
n = ceil(carlen/2);
newvec = conv([oldvec(1:n) oldvec oldvec(end-n+1:end)], car, 'same');
newvec(1:n) = []; newvec(end-n+1:end) = [];
