function [xQ, xE, codebook] = xquant(x, nbits, maxValue, method)
% least significant bit for nominal range [-1 ... 1]
lsb = 2^-(nbits-1)*maxValue/2;

% ---------------------- quantize signal ---------------------- %
    switch method

        % mid-tread quantisation, with a value for zero but one quantisation
        % value less for values >07

        case 'mid-tread'
            % partition defines boundaries used for quantisation
            partition = -1+.5*lsb:lsb:1-1.5*lsb;
            % codebook defines values used for quantisation
            codebook  = -1:lsb:1-lsb;

        % mid-rise quantisation. No value for zero but evan spacing of values
        % around zero
        case 'mid-rise'
            partition = lsb:lsb:maxValue-lsb;
            codebook  = 0.5*lsb:lsb:maxValue-.5*lsb;
    end
[~, xQ] = quantiz(x,partition,codebook);
% calculate quantisation error

%xE = xQ - x;
xQ = xQ';
%xE = xE';
end