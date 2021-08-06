function spgdemo_seismic_irr()
% fft work with inverse interpolation
% curvelet work with matrix transpose
% stft does not converge
    addpath(genpath('../spgl1-2.1'));
    addpath(genpath('../barylag1d'));
    
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    % Create random m-by-n encoding matrix and sparse vector
    n = 256; n2 = 128; m = 96;
%   n = 512;  n2 = 512; m = 256;
    rp = randperm(n2);
    k = 3;
    load ('../barylag1d/seismic.mat');
    f0 = (data{2}(1:n,1:n2)-128)/128;
    
    l = (1:n2)';
    p = l(sort(rp(1:m),'ascend'));
    
%     p = sort(rand(m,1)*n2,'ascend');
    
    p(end) = l(end); p(1) = l(1);
    
    dh = 0.05;
    p = p + dh;

    [A, AT] = construct_sampler(k,l,p);

    % Set up vector b, and run solver
    b = A (f0) + randn(n,m) * 0.00;
    sigma = norm(f0)*0.02;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);
    
    v2m_b = @(x) reshape(x,[n,m]);
    m2v_b = @(x) x(:);
    
%     [Psi,PsiT,v2m_x,m2v_x] = construct_FFT(n,n2);    % FFT    
%     [Psi,PsiT,v2m_x,m2v_x] = construct_STFT(n,n2);  % STFT
    [Psi,PsiT,v2m_x,m2v_x] = construct_Curvelet(n,n2);  % Curvelet
    opts.weights = m2v_x(construct_Curvelet_weight(n,n2));
    
    Ah = @ (x,mode)  Measure(A,AT,Psi,PsiT,x,v2m_x,v2m_b,m2v_x,m2v_b,mode);

    xv = spg_bpdn(Ah, b(:), sigma, opts);
    
    f = PsiT(v2m_x(xv));
    
    c = max(abs(f0(:)));
    subplot(121);imagesc(b,[-c,c]);colormap gray;
    subplot(122);imagesc(real(f),[-c,c]);colormap gray;
    
end

function y = Measure(A,AT,Psi,PsiT,x,v2m_x,v2m_b,m2v_x,m2v_b,mode)
    if mode==1
        y = A(PsiT(v2m_x(x)));
        y = m2v_b(y);
    elseif mode==2
        y = Psi(AT(v2m_b(x)));
        y = m2v_x(y);
    end
end

function [A, AT] = construct_sampler(k,l,p)
%     A = @(x) barylag_k_vec(k,l,x,p); 
%     AT= @(x) barylag_k_vec(k,p,x,l); % inverse interpolation
    A = @(x) barylag_k_mat2d(k,l,x,p,1);
    AT= @(x) barylag_k_mat2d(k,l,x,p,2); % matrix transpose 
%     AT= @(x) barylag_k_mat2d(k,p,x,l,1);   % inverse interpolation
end
    
function [Psi,PsiT,v2m_x,m2v_x] = construct_STFT(nt,nx)
    addpath(genpath('../gwt'));    
    noverlap = 15;
    nwindow  = 192;
    dt = 0.004;
    dx = 20;
    Psi = @(x) ystft2d(x,noverlap,nwindow,dt,dx);
    PsiT= @(x) iystft2d(x,noverlap,nwindow,dt,nt,nx);
    s = size(Psi(ones(nt,nx)));
    v2m_x = @(x) reshape(x,s);
    m2v_x = @(x) x(:);
end

function [Psi,PsiT,v2m_x,m2v_x] = construct_FFT(nt,nx)
    Psi = @(x) fft2(x)/sqrt(numel(x));
    PsiT = @(x)ifft2(x)*sqrt(numel(x));
    v2m_x = @(x) reshape(x,[nt,nx]);
    m2v_x = @(x) x(:);
end

function [Psi,PsiT,v2m_x,m2v_x] = construct_Curvelet(nt,nx)
    addpath(genpath('../CurveLab-2.1.3/fdct_wrapping_matlab'));
    Psi = @(x) fdct_wrapping(x,1,2);
    PsiT = @(x) ifdct_wrapping(x,1);
    
    d = zeros(nt,nx);
    c = Psi(d);
    nc = cellnum(c);
    m2v_x = @(x) cell2vec(x,nc);
    v2m_x = @(x) vec2cell(x,c);
end

function E = construct_Curvelet_weight(nt,nx)
    F = ones(nt,nx);
    X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
    C = fdct_wrapping(X,0,2);
    % Compute norm of curvelets (exact)
    E = cell(size(C));
    for s=1:length(C)
      E{s} = cell(size(C{s}));
      for w=1:length(C{s})
        A = C{s}{w};
        E{s}{w} = ones(size(A))*sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
      end
    end

end

function C=vec2cell(v,C)
    k = 1;
    for s=1:length(C)
      for w=1:length(C{s})
          n1 = numel(C{s}{w});
            C{s}{w}(:) = v(k:k+n1-1);
          k = k+n1;
      end
    end
end

function nc = cellnum(C)
    nc = 0;
    for s=1:length(C)
      for w=1:length(C{s})
          nc = nc+ numel(C{s}{w});
      end
    end
end

function cfs = cell2vec(C,nc)
    cfs =zeros(nc,1);
    k = 1;
    for s=1:length(C)
      for w=1:length(C{s})
          n1 = numel(C{s}{w});
          cfs(k:k+n1-1) = C{s}{w}(:);
          k = k+n1;
      end
    end
end
