function spgdemo_dong_sea()
    addpath(genpath('../spgl1-2.1'));
    addpath(genpath('../barylag1d'));
    addpath(genpath('../SeismicLab'));
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    m = 60;
    k = 3;
%     A(:,1:2:end) = 1;
    load ('sea.mat');
    f0 = gain(d,0.004,'agc',0.5,1);
    f0(isnan(f0)) = 0;
    imagesc(f0);colormap gray
    [n,n2] = size(f0);
    
    l = (1:n2)';
    
    rp = randperm(n2);
    p = l(sort(rp(1:m),'ascend'));
    
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
    Ah = @ (x,mode)  Measure(A,AT,x,v2m_x,v2m_b,m2v_x,m2v_b,mode);

    xv = spg_bpdn(Ah, m2v_b(b), sigma, opts);
    
    f = ifft2(v2m_x(xv))*sqrt(numel(f0));

    subplot(121);imagesc(b,[-1,1]);colormap gray;
    subplot(122);imagesc(real(f),[-1,1]);colormap gray;
    
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


