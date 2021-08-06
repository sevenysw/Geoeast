function spgdemo_image()
    addpath(genpath('../spgl1-2.1'));
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    % Create random m-by-n encoding matrix and sparse vector
    n2 = 10;
    m = 50; n = 128; k = 14*n2;
    [A,~] = qr(randn(n,m),0);
    A  = A';
    
    p  = randperm(n*n2); p = p(1:k);
    x0 = zeros(n,n2); x0(p) = randn(k,1);

    % Set up vector b, and run solver
    b = A * x0 + randn(m,n2) * 0.075;
    sigma = 0.10;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);

    % A as a function handle
    m2v_x = @(x) reshape(x,[n*n2,1]);
    v2m_x = @(x) reshape(x,[n,n2]);
    
    m2v_b = @(x) reshape(x,[m*n2,1]);
    v2m_b = @(x) reshape(x,[m,n2]);
    Ah = @ (x,mode)  Measure(A,x,v2m_x,v2m_b,m2v_x,m2v_b,mode);

    xv = spg_bpdn(Ah, m2v_b(b), sigma, opts);
    
    x = v2m_x(xv);

    subplot(121);imagesc(x0);
    subplot(122);imagesc(x);
end

function y = Measure(A,x,v2m_x,v2m_b,m2v_x,m2v_b,mode)
    if mode==1
        y = A*v2m_x(x);
        y = m2v_b(y);
    elseif mode==2
        y = A'*v2m_b(x);
        y = m2v_x(y);
    end
end


