f = @(x) x.^2;
tic
x = randn(10000,1);
nt = 10000;
for i=1:nt
    y = f(x);
end
toc

tic
nt = 10000;
for i=1:nt
    f = @(x) x.^2;
    y = f(x);
end
toc