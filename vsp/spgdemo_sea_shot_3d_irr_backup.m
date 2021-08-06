function spgdemo_seismic_irr()
    addpath(genpath('../spgl1-2.1'));
    addpath(genpath('../barylag1d'));
    addpath(genpath('../commonfunction'));
    addpath(genpath('../SeismicLab'));

    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    k = 3;
    load ('sea3d_2.mat');
    f0 = D(257:512,:,1:128);
    fshow = f0*0;
    [nt,n2,n3] = size(f0);
    
%     for i=1:n3
%         f0(:,:,i) = gain(f0(:,:,i),0.004,'agc',0.5,1);
%     end
    
    f0(isnan(f0)) = 0;
        
    l = x2'+0.0005; % original, np
    p = x1'; % sampled, n2
    np = length(x2);
%     A = @(x) barylag_k_mat(k,l,x,p);
%     AT= @(x) barylag_k_mat(k,p,x,l);

    A = @(x) barylag_k_mat3d(k,l,x,p,1);
    AT= @(x) barylag_k_mat3d(k,p,x,l,1);
    
    b = f0;
    
    sigma = norm(f0(:))*0.020;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);

    % A as a function handle
    v2m_x = @(x) reshape(x,[nt,np,n3]);
    
    v2m_b = @(x) reshape(x,[nt,n2,n3]);
    
    Ah = @ (x,mode)  Measure(A,AT,x,v2m_x,v2m_b,mode);

    xv = spg_bpdn(Ah, b(:), sigma, opts);
    
    f = real(ifftn(v2m_x(xv))*sqrt(numel(f0)));
    
    c = 100;
    seishow3D(f0,100,-c,c);
    seishow3D(f,100,-c,c);
    figure,imagesc(reshape(permute(f0,[1 3 2]),[nt,n2*n3]),[-c,c]);colormap gray;
    figure,imagesc(reshape(permute(f,[1 3 2]),[nt,np*n3]),[-c,c]);colormap gray;
end

function y = Measure(A,AT,x,v2m_x,v2m_b,mode)
    if mode==1
        y = A(ifftn(v2m_x(x))*sqrt(numel(x)));
    elseif mode==2
        x = AT(v2m_b(x));
        y = fftn(x)/sqrt(numel(x)); 
    end
    y = y(:);
end

