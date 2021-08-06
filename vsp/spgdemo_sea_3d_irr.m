function spgdemo_seismic_irr()
    addpath(genpath('../spgl1-2.1'));
    addpath(genpath('../barylag1d'));
    addpath(genpath('../commonfunction'));

    % Initialize random number generators 
    rand('state',0);
    randn('state',0);

    k = 3;
    load ('sea3d.mat');
    f0 = D(128:512,:,:);
    fshow = f0*0;
    [nt,n2,n3] = size(f0);
    
    for i=1:n3
        f0(:,:,i) = gain(f0(:,:,i),0.004,'agc',0.5,1);
    end
    
    f0(isnan(f0)) = 0;
    
    xs = [nt,n2,n3];
    
    l = (1:n2)';
    m0 = 80;
    
    p = {};
    np = [];
    for i=1:n3     
%         p{i} = sort(rp(1:round(n2/2+randn*2)),'ascend');

        rp = randperm(n2)';
        m = round(m0+randn*2);
        p{i} = l(sort(rp(1:m),'ascend'));
        
        p{i}(end) = l(end); p{i}(1) = l(1);
        np(i) = length(p{i});
    end

    for i=1:n3
        A{i} = @(x) barylag_k_vec(k,l,x,p{i});
        AT{i}= @(x) barylag_k_vec(k,p{i},x,l);
        b{i} = A{i}(f0(:,:,i));
        fshow(:,1:np(i),i) = b{i};
    end
    
    sigma = 0.020;       % Desired ||Ax - b||_2
    opts = spgSetParms('verbosity',1);

    % A as a function handle
    v2m_x = @(x) reshape(x,xs);
    
    m2v_b = @(x) cell2vec(x);
    v2m_b = @(x) vec2cell(x,nt,np);
    
    Ah = @ (x,mode)  Measure(A,AT,x,v2m_x,v2m_b,m2v_b,xs,mode);

    xv = spg_bpdn(Ah, m2v_b(b), sigma, opts);
    
    f = real(ifftn(v2m_x(xv))*sqrt(numel(f0)));
    
    c = 1;
    seishow3D(f0,100,-c,c);
    seishow3D(fshow,100,-c,c);
    seishow3D(f,100,-c,c);
end

function y = Measure(A,AT,x,v2m_x,v2m_b,m2v_b,xs,mode)
    if mode==1
        x = ifftn(v2m_x(x))*sqrt(numel(x));
        for i=1:xs(3)
            y{i} = A{i}(x(:,:,i));
        end
        y = m2v_b(y);
    elseif mode==2
        b = v2m_b(x);
        x = zeros(xs);
        for i=1:length(b)    
            x(:,:,i) = AT{i}(b{i});
        end
        y = fftn(x)/sqrt(numel(x));
        y = y(:);
    end
end

function y = cell2vec(x)
    y = [];
    for i=1:length(x)
        y = [y;x{i}(:)];
    end
end

function y = vec2cell(x,nt,np)
    y = {};
    s = 1;
    for i=1:length(np)
        mat = reshape(x(s:s+nt*np(i)-1),[nt,np(i)]);
        y{i} = mat;
        s = s+nt*np(i);
    end
end
