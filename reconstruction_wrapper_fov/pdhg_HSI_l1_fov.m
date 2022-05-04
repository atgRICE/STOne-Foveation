function [x,outs ] = pdhg_HSI_l1_fov( R, b, mu, specWeight, order, res_prev, topLeft, botRight, record_ds, fovParams, guess )

order_prev = createOrderingData(res_prev,'full');
numBands = size(b,2);

%%  Allocate memory
rows = order.n1;
cols = order.n2;
x0=zeros(cols,numBands);
y0=zeros(3*cols+rows,numBands);


%% Get Derivative Matrices
% Get linear indices of foveated region
[Dx,Dy] = createDifferenceOperators_fov(order,order_prev,topLeft,botRight);
Dxt = Dx';
Dyt = Dy';

tic
scale = max(max(abs(b)));
coeff = 100;
b = b*coeff/scale;

xProx =@(x,tau) x;
yProx = @(y,sigma)[projectInf(y(1:cols,:),y(cols+1:2*cols,:) ) ;...
                   min(max(y(2*cols+1:3*cols,:),-1),1); ...
                   min(max(y(3*cols+1:end,:)-sigma*(R.*b),-mu),mu) ];
w = [1,1,specWeight];
A = @(x) [w(1)*Dx*x; 
          w(2)*Dy*x; 
          w(3)*DDz(x); 
          R.*STO(UndoRLE(x,record_ds,fovParams))];
At = @(y) w(1)*Dxt*y(1:cols,:)+...
          w(2)*Dyt*y(cols+1:2*cols,:)+...
          w(3)*DDzt(y(2*cols+1:3*cols,:))+...
          fSTO_RLE(R.*y(3*cols+1:end,:),fovParams) ;

%%  Determine timestep paramters 
opts = [];
%opts.L = 0.95/8;
opts.tol = 0.1;
opts.maxIters = 3000;
opts.verbose = false; 

if exist('guess','var')
x0 = guess;
y0 = [sign(Dx*x0) ; sign(Dy*x0); sign(R.*fSTO_RLE(UndoRLE(x0,record_ds,rows/res_prev^2),fovParams)')];
end

tic
[x,outs]= pdhg_adaptive(x0,y0,A,At,xProx, yProx,opts);            
% outs.times = toc;
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