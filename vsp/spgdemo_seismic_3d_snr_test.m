% whole data SNR test
function spgdemo_seismic()
    addpath(genpath('../spgl1-2.1'));
    addpath(genpath('../commonfunction'));
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    % Create random m-by-n encoding matrix and sparse vector

    load ('/Volumes/Yu/BaiduYun/dong_data/GulfOfSuez128.mat');
    
    f0 = D;
    [nt,n2,n3] = size(f0);
    M2 = zeros(n2,n3);
    rp = randperm(n2*n3);
    M2(rp(1:n2*n3/4)) = 1;
    M = repmat(reshape(M2,[1,n2,n3]),[nt,1,1]);
%     A(:,1:2:end) = 1;    
    % Set up vector b, and run solver
    A = @(x)M.*x;
    AT = @(x)M.*x;
    b = A (f0) + randn(nt,n2,n3) * 0.00;
    sigma = 0.020;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);

    % A as a function handle
    v2m_x = @(x) reshape(x,[nt,n2,n3]);
    v2m_b = @(x) reshape(x,[nt,n2,n3]);
    
    Ah = @ (x,mode)  Measure(A,AT,x,v2m_x,v2m_b,mode);

    xv = spg_bpdn(Ah, b(:), sigma, opts);
    
    f = real(ifftn(v2m_x(xv))*sqrt(numel(f0)));

    seishow3D(b);
    seishow3D(f);
    
    ysnr = @(x,y) 20 * log(norm(x(:))/norm((x(:)-y(:))))/log(10);
    ysnr(f,f0)
    
end

function y = Measure(A,AT,x,v2m_x,v2m_b,mode)
    if mode==1
        y = A(ifftn(v2m_x(x)))*sqrt(numel(x));
    elseif mode==2
        y = fftn(AT(v2m_b(x)))/sqrt(numel(x));
    end
    y = y(:);
end


