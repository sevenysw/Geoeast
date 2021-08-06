function spgdemo_image_fft()
    addpath(genpath('../spgl1-2.1'));
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    % Create random m-by-n encoding matrix and sparse vector
    n2 = 100;
    m = 50; n = 128; k = 4;
    [A,~] = qr(randn(n,m),0);
    A  = A';
%     A  = randn(m,n);
    load ../barylag1d/seismic.mat
    p  = randperm(n*n2); p = p(1:k);
    x0 = zeros(n,n2); x0(p) = randn(k,1);
    f0 = (ifft2(x0))*sqrt(numel(x0));
    
    f0 = data{1}(1:n,1:n2);

    % Set up vector b, and run solver
    b = A * f0 + randn(m,n2) * 0.00;
    sigma = 0.20;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);

    % A as a function handle
    m2v_x = @(x) reshape(x,[n*n2,1]);
    v2m_x = @(x) reshape(x,[n,n2]);
    
    m2v_b = @(x) reshape(x,[m*n2,1]);
    v2m_b = @(x) reshape(x,[m,n2]);
    Ah = @ (x,mode)  Measure(A,x,v2m_x,v2m_b,m2v_x,m2v_b,mode);

    xv = spg_bpdn(Ah, m2v_b(b), sigma, opts);
    
    f = ifft2(v2m_x(xv))*sqrt(numel(x0));

    subplot(121);imagesc(abs(f0));
    subplot(122);imagesc(abs(f));
    
end

function y = Measure(A,x,v2m_x,v2m_b,m2v_x,m2v_b,mode)
    if mode==1
        y = A*ifft2(v2m_x(x))*sqrt(numel(x));
        y = m2v_b(y);
    elseif mode==2
        y = fft2(A'*v2m_b(x))/sqrt(numel(x));
        y = m2v_x(y);
    end
end


