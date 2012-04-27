%{
    Can BP solve a LP in standard form?
    Yes, if we smooth

(P)     min c'x s.t. x >= 0, Ax=b

when we form the Lagrangian, since everythingn is linear, when we minimize
over x, we get either 0 or -Inf.  Adding a quadratic term changes this!
And that means the two dual variables, lambda and nu, are no longer coupled.

Traditional dual:

(D)    max -b'*nu
    s.t. lambda >= 0, and lambda = c + A'*nu

If we smooth c'x to c'x + mu/2*|||x-x0||^2, dual is now uncoupled quadratic


%}

%{
N = 3500;
M = round(N/2);
x = randn(N,1)+10;
% A = randn(M,N);
A = sprand(M,N,.05);
c = randn(N,1);
b = A*x;

%}

cd ~/Documents/research/software/cvx;
cvx_setup;
cd ~/Documents/research/software/scallop;
load LP.A;
A = spconvert(LP);
clearvars LP;
b = importdata('LP.b');
c = importdata('LP.c');
N = length(c);

cvx_begin
    cvx_precision low
    variable x(N)
    minimize( c'*x )
    subject to
        x >= 0
        A*x == b
cvx_end
x_cvx = x;
c'*x_cvx
cvx_status
exit;



%% solve in TFOCS

% obj    = smooth_linear(c)
%{
opts = [];
opts.maxIts     = 4000;
opts.errFcn        = {@(f,d,x) norm( x - x_cvx )/norm(x_cvx )};
opts.errFcn{end+1} = @(f,d,x) norm( A*x - b )/norm(b );
opts.restart    = 1500;
opts.continuation = true;
% x0   = [];
%x0 = A\b;
x0 = zeros(size(A,2),1);
x0 = sign(x0).*max( abs(x0) - norm(x0,Inf)/10, 0 );
z0   = [];
mu = 1e-1;

opts.stopCrit   = 4;
opts.tol        = 1e-8;
contOpts        = [];       % options for "continuation" settings
%   (to see possible options, type "continuation()" )
contOpts.maxIts = 5;
contOpts.accel  = false;
% contOpts.betaTol = 1.3; % how much to decrease tolerance
contOpts.betaTol    = 10;
% contOpts.innerMaxIts = 400;
contOpts.finalTol = 1e-10;

[x,out,optsOut] = solver_sLP(c,A,b, mu, x0, z0, opts, contOpts);

%% plot
%figure(2);
%semilogy( out.err(:,1) );
%hold all
%emphasizeRecent
%}
