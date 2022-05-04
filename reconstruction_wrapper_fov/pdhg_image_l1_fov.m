function [x,outs ] = pdhg_image_l1_fov( R, b, mu, w, order, res_prev, topLeft, botRight, record_ds, fovParams, guess )

order_prev = createOrderingData(res_prev,'full');

%%  Allocate memory
rows = order.n1;
cols = order.n2;
x0=zeros(cols,1);
y0=zeros(2*cols+rows,1);

%% Get Derivative Matrices
[Dx,Dy] = createDifferenceOperators_fov(order,order_prev,topLeft,botRight);
Dxt = Dx';
Dyt = Dy';

scale = max(max(abs(b)));
b = b*100/scale;

xProx =@(x,tau) x;
yProx = @(y,sigma)[projectInf(y(1:cols,:),y(cols+1:2*cols,:) ) ; min(max(y(2*cols+1:end,:)-sigma*(R.*b),-mu),mu) ];
A = @(x) [w(1)*Dx*x ; 
          w(2)*Dy*x; 
          R.*STO(UndoRLE(x,record_ds,fovParams))];
At = @(y) w(1)*Dxt*y(1:cols,:)+...
          w(2)*Dyt*y(cols+1:2*cols,:)+...
          fSTO_RLE(R.*y(2*cols+1:end,:),fovParams);
      
%%  Determine timestep paramters 
opts = [];
%opts.L = 0.95/8;
opts.tol = .01;
opts.maxIters = 3000;
opts.verbose = false; 

if exist('guess','var')
x0 = guess;
y0 = [sign(Dx*x0) ; sign(Dy*x0); sign(R.*fSTO_RLE(UndoRLE(x,record_ds,rows/res_prev^2),fovParams)')];
end

[x,outs]= pdhg_adaptive(x0,y0,A,At,xProx, yProx,opts);
x = x*scale/100;
