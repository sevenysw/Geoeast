function spgdemo_seismic()
    addpath(genpath('../spgl1-2.1'));
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    % Create random m-by-n encoding matrix and sparse vector
    n2 = 128;
    n = 256; 
    M = zeros(n,n2);
    rp = randperm(n2);
    M(:,rp(1:64)) = 1;
%     A(:,1:2:end) = 1;
    load ('../barylag1d/seismic.mat');
    
    f0 = data{4}(1:n,1:n2)-128;
    % Set up vector b, and run solver
    A = @(x)M.*x;
    AT = @(x)M.*x;
    b = A (f0) + randn(n,n2) * 0.00;
    sigma = 0.020;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);

    % A as a function handle
    m2v_x = @(x) reshape(x,[n*n2,1]);
    v2m_x = @(x) reshape(x,[n,n2]);
    
    m2v_b = @(x) reshape(x,[n*n2,1]);
    v2m_b = @(x) reshape(x,[n,n2]);
    Ah = @ (x,mode)  Measure(A,AT,x,v2m_x,v2m_b,m2v_x,m2v_b,mode);

    xv = spg_bpdn(Ah, m2v_b(b), sigma, opts);
    
    f = ifft2(v2m_x(xv))*sqrt(numel(f0));

    subplot(121);imagesc(b,[-128,128]);colormap gray;
    subplot(122);imagesc(real(f),[-128,128]);colormap gray;
    
end

function y = Measure(A,AT,x,v2m_x,v2m_b,m2v_x,m2v_b,mode)
    if mode==1
        y = A(ifft2(v2m_x(x)))*sqrt(numel(x));
        y = m2v_b(y);
    elseif mode==2
        y = fft2(AT(v2m_b(x)))/sqrt(numel(x));
        y = m2v_x(y);
    end
end


