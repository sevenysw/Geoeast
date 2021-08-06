function spgdemo_seismic_irr()
    addpath(genpath('../spgl1-2.1'));
    addpath(genpath('../barylag1d'));
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    % Create random m-by-n encoding matrix and sparse vector
    n2 = 128;
    n = 256; 
    m = 64;
    rp = randperm(n2);
    k = 3;
%     A(:,1:2:end) = 1;
    load ('../barylag1d/seismic.mat');
    f0 = (data{2}(1:n,1:n2)-128)/128;
    
    l = (1:n2)';
%     p = l(sort(rp(1:m),'ascend'));
    
    p = sort(rand(m,1)*n2,'ascend');
    
    p(end) = l(end); p(1) = l(1);
    
    A = @(x) barylag_k_vec(k,l,x,p);
    AT= @(x) barylag_k_vec(k,p,x,l);
    
    % Set up vector b, and run solver
    b = A (f0) + randn(n,m) * 0.00;
    sigma = 0.020;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);

    % A as a function handle
    m2v_x = @(x) reshape(x,[n*n2,1]);
    v2m_x = @(x) reshape(x,[n,n2]);
    
    m2v_b = @(x) reshape(x,[n*m,1]);
    v2m_b = @(x) reshape(x,[n,m]);
    xsize = [n,n2];
    bsize = [n,m];
 
    Ah = @ (x,mode)  Wavedb1(A,AT,x,xsize,bsize,mode);

    xv = spg_bpdn(Ah, m2v_b(b), sigma, opts);
    
    level = 3;wname = 'db1';
    [~,s] = wavedec2(reshape(xv,xsize),level,wname);
    f = waverec2(xv,s,wname);

    subplot(121);imagesc(b,[-1,1]);colormap gray;
    subplot(122);imagesc(real(f),[-1,1]);colormap gray;
    
end


function y = Wavedb1(A,AT,x,xsize,bsize,trans)
persistent s;
level = 3;
wname = 'db1';

if trans==1 
    [~,s] = wavedec2(reshape(x,xsize),level,wname);
    Y = A(waverec2(x,s,wname));
else       
    [Y,s] = (wavedec2(AT(reshape(x,bsize)),level,wname));
end
y = Y(:);
end



function y = Measure(A,AT,x,v2m_x,v2m_b,m2v_x,m2v_b,mode)
    if mode==1
        y = A(ifft2(v2m_x(x))*sqrt(numel(x)));
        y = m2v_b(y);
    elseif mode==2
        y = fft2(AT(v2m_b(x)))/sqrt(numel(x));
        y = m2v_x(y);
    end
end
