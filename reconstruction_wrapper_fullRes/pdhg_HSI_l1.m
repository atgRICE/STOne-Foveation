function [x,outs ] = pdhg_HSI_l1( R, b, mu, SpecWeight, order, guess )

rows = order.n; 
cols = order.n;
numBands = size(R,2);

%%  Allocate memory
NN = rows*cols;
x0=zeros(NN,numBands);
y0=zeros(4*NN,numBands);


%% Get Derivative Matrices
[Dx,Dy] = createDifferenceOperators(order);
Dxt = Dx';
Dyt = Dy';

scale = max(max(abs(b)));
coeff = 100;
b = coeff*b/scale;

xProx =@(x,tau) x;
yProx = @(y,sigma)[projectInf(y(1:end/4,:),y(end/4+1:2*end/4,:) ) ; min(max(y(2*end/4+1:3*end/4,:),-1),1) ; min(max(y(3*end/4+1:end,:)-sigma*(R.*b),-mu),mu) ];

weights = [1 1 SpecWeight];
A = @(x) [weights(1)*Dx*x ; weights(2)*Dy*x; weights(3)*DDz(x) ; R.*STO(x)];
At = @(y) weights(1)*Dxt*y(1:end/4,:)+weights(2)*Dyt*y(end/4+1:2*end/4,:)+weights(3)*DDzt(y(2*end/4+1:3*end/4,:))+STO(R.*y(3*end/4+1:end,:)) ;

%%  Determine timestep paramters 
opts = [];
opts.tol = .1;
opts.maxIters = 2400;
opts.verbose = false; 

if exist('guess','var')
x0 = guess;
y0 = [sign(Dx*x0) ; sign(Dy*x0); sign(Dz*x0) ; sign(R.*STO(x0))];
end

tic
[x,outs]= pdhg_adaptive(x0,y0,A,At,xProx, yProx,opts);

x = x*scale/coeff;

end

function dz = Dz(u)
dz = imfilter(u, [-1 1 0]);
dz(:,1) = 0;
end

function dz = Dzt(u)
dz = imfilter(u, [0 1 -1]);
dz(:,1) = -u(:,2);
dz(:,end) = u(:,end);
end

function dz = DDz(u)
dz = Dz(Dz(u));
end

function dz = DDzt(u)
dz = Dzt(Dzt(u));
end