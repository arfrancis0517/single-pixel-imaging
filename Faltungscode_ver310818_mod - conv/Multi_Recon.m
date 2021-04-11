
figure
tic
numFrames=15000;
smasks=(numFrames/2)-1;

period1=15000;
skip1=30000;
skip2=0;

data2(1)=data(1);
data2(2)=data(2);

    for i=3:period1;    
        data2(i)=data(i+skip1);
    end

    for i=period1+1:numFrames;    
        data2(i)=data(i+skip2);
    end



    for i=1:period1/2-1;    
        indices2(i)=indices(i+skip1/2);
    end
    for i=period1/2:smasks;    
        indices2(i)=indices(i+skip2/2);
    end



%%

%Read data from white masks

s2=3;
for i=1:(numFrames-smasks)-2;
    aux_w(i)=data2(s2);
    s2=s2+2;
end

%%

%Read data from structured masks
ym=zeros(1);
s4=4;
for i=1:smasks;      
    ym(i)=data2(s4);
    s4=s4+2;
end
%----------------

%
% Driftcorrection

for i=1:smasks;
ym(i)=ym(i).*aux_w(400)./aux_w(i);
end
for i=1:smasks;
aux_w(i)=aux_w(400);
end

%----------------
%}

%%
% Reconstruction options

opt = struct;
opt.tol = 1e-6; %tolerance
opt.max_time = 10 * 60; %maximum calculation time in seconds
opt.max_iter = 10000; %maximimum iterations


opt.method = 'flexbox';
opt.flexbox_internal_iter = 100;  %How many iterations before checkpoint -> callback
opt.flexbox_nnls = 1; %use nnls or not
opt.flexbox_l1grad = 1e3; %factor for ||grad(f)||_1 regularization (this is TV-1-Norm)
opt.flexbox_l2grad = 0.00; %factor for ||grad(f)||_2^2 regularization (TV-2-Norm)
opt.flexbox_wname = 'none'; %use some wavelet or 'none'

[x, mu] = reconstruct(kernel, indices2', ym', opt);
imagesc(x);
toc
