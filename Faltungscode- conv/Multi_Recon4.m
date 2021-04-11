
figure
tic
numFrames=15000;
n=40;    %number of uploaded streams
smasks=(numFrames)*3/4;


%%

%For multiply upload streams
for i=1:(n-1);    
    data(15000*i+3)=data(15000*i+7);
end
%%

add=zeros(1080,1920);

for m=0:n-1;
    
period1=15002;
skip1=15000*m;
skip2=0;

data2(1)=data(1);
data2(2)=data(2);

for i=3:period1;    
        data2(i)=data(i+skip1);
end
    
for i=1:smasks;    
        indices2(i)=indices(i+11250*m);
    end
    
%{
    for i=2:period1;    
        data2(i)=data(i+skip1+2);
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
%}


%%

%Read data from white masks

s2=3;
state=0;
for i=1:smasks;
     
    if (state==4)
        s2=s2+4;
        aux_w(i)=data2(s2);
        state=1;        
     else
        aux_w(i)=data2(s2);
        state=state+1;        
    end
     
end

%%

%Read data from structured masks
ym=zeros(1);
s4=3;
state=1;
for i=1:smasks;      
    
    
    if (state==4)
        s4=s4+2;
        ym(i)=data2(s4);
        state=2;        
    else
        s4=s4+1;
        ym(i)=data2(s4);        
        state=state+1;        
    end
    
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
add=add+x;
%imagesc(x)
%figure
end
imagesc(add);

toc
