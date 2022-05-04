function [x,outs ] = pdhg_HSVid_l1_fovea_squareDeriv_stable_multiROIsingleMat( R, b, mu, numBands, numFrames, order, res_prev, topLeft, botRight, record_ds, fovParams, guess )

order_prev = createOrderingData(res_prev,'full');

%%  Allocate memory
rows = order.n1;
cols = order.n2;
x0=zeros(cols,numBands*numFrames);
y0=zeros(4*cols+rows,numBands*numFrames);

%% Get Derivative Matrices
% Get linear indices of foveated region
[Dx,Dy] = createDifferenceOperators_fov(order,order_prev,topLeft,botRight);
Dxt = Dx';
Dyt = Dy';

scale = max(max(abs(b)));
b = b*100/scale;

xProx =@(x,tau) x;
yProx = @(y,sigma) [projectInf(y(1:cols,:),y(cols+1:2*cols,:) ) ;...
                    min(max(y(2*cols+1:3*cols,:),-1),1); ...
                    min(max(y(3*cols+1:4*cols,:),-1),1); ...
                    min(max(y(4*cols+1:end,:)-sigma*(R.*b),-mu),mu) ];
w = [0.7, 0.7, 0.5, 0.3];
A = @(x) [w(1)*Dx*x ; w(2)*Dy*x; ...
          w(3)*DDz_new(x,numBands,numFrames); ...
          w(4)*Dt_new(x,numBands,numFrames); ...
          R.*STO(UndoRLE(x,record_ds,fovParams))];
      
At = @(y) w(1)*Dxt*y(1:cols,:)+w(2)*Dyt*y(cols+1:2*cols,:)+...
          w(3)*DDzt_new(y(2*cols+1:3*cols,:),numBands,numFrames)+...
          w(4)*Dtt_new(y(3*cols+1:4*cols,:),numBands,numFrames)+...
          fSTO_RLE(R.*y(4*cols+1:end,:),fovParams);
      
%%  Determine timestep paramters 
opts = [];
%opts.L = 0.95/8;
opts.tol = .01;
opts.maxIters = 3000;
opts.verbose = true; 

if exist('guess','var')
x0 = guess;
y0 = [sign(Dx*x0) ; sign(Dy*x0); sign(R.*STOf_poc(x0,bigpat))];
end

tic
[x,outs]= pdhg_adaptive_nonneg(x0,y0,A,At,xProx, yProx,opts);            
% outs.times = toc;
x = x*scale/100;

end

function dz = Dz_new(u, numBands, numFrames)
dz = zeros(size(u));
for ii = 1:numFrames
    dz(:,((ii-1)*numBands+1):ii*numBands) = Dz(u(:,((ii-1)*numBands+1):ii*numBands));
end
end

function dz = Dz(u)
dz = imfilter(u, [-1 1 0]);
dz(:,1) = 0;
end

function dt = Dt_new(u, numBands, ~)
dt = zeros(size(u));
for ii = 1:numBands
    dt(:,ii:numBands:end) = Dz(u(:,ii:numBands:end));
end
end

function dz = Dzt_new(u, numBands, numFrames)
dz = zeros(size(u));
for ii = 1:numFrames
    dz(:,((ii-1)*numBands+1):ii*numBands) = Dzt(u(:,((ii-1)*numBands+1):ii*numBands));
end
end

function dz = Dzt(u)
dz = imfilter(u, [0 1 -1]);
dz(:,1) = -u(:,2);
dz(:,end) = u(:,end);
end

function dt = Dtt_new(u, numBands, ~)
dt = zeros(size(u));
for ii = 1:numBands
    dt(:,ii:numBands:end) = Dzt(u(:,ii:numBands:end));
end
end

function ddz = DDz_new(u, numBands, numFrames)
ddz=Dz_new(Dz_new(u, numBands, numFrames), numBands, numFrames);
end

function ddz = DDzt_new(u, numBands, numFrames)
ddz=Dzt_new(Dzt_new(u, numBands, numFrames), numBands, numFrames);
end
