function spgdemo()
    addpath(genpath('../spgl1-2.1'));
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    % Create random m-by-n encoding matrix and sparse vector
    m = 50; n = 128; k = 14;
    [A,Rtmp] = qr(randn(n,m),0);
    A  = A';
    p  = randperm(n); p = p(1:k);
    x0 = zeros(n,1); x0(p) = randn(k,1);

    % -----------------------------------------------------------
    % Solve the basis pursuit denoise (BPDN) problem:
    %
    %    minimize ||x||_1 subject to ||Ax - b||_2 <= 0.1
    %
    % -----------------------------------------------------------

    % Set up vector b, and run solver
    b = A * x0 + randn(m,1) * 0.075;
    sigma = 0.10;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);

    % A as a function handle
    Ah = @ (x,mode)  Measure(A,x,mode);

    x = spg_bpdn(Ah, b, sigma, opts);

    figure(1);
    plot(1:n,x,'b', 1:n,x0,'ro');
    legend('Recovered coefficients','Original coefficients');
    title('(b) Basis Pursuit Denoise');
end

function y = Measure(A,x,mode)
    if mode==1
        y = A*x;
    elseif mode==2
        y = A'*x;
    end
end


