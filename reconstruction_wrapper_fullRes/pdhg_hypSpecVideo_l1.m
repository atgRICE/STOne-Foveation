function [x,outs] = pdhg_hypSpecVideo_l1( R, b, mu, numBands, numFrames, order )

rows = order.n; 
cols = order.n;

%%  Allocate memory
NN = rows*cols;
x0 = zeros(NN,numBands*numFrames);
y0 = zeros(5*NN,numBands*numFrames);  %  x,y,z,t store the derivatives in the 4 dimensions


%% Get Derivative Matrices
[Dx Dy] = createDifferenceOperators(order);
Dxt = Dx';
Dyt = Dy';

%b = b*10000/max(max(abs(b))); % 10000
scale = max(b(:))-min(min(b(b~=0)));
% b = b-mean(b(:));
b = 100*b/scale;

xProx =@(x,tau) x;
yProx = @(y,sigma) [projectInf(y(1:end/5,:),y(end/5+1:2*end/5,:)) ; ...
                    min(max(y(2*end/5+1:3*end/5,:),-1),1); ...
                    min(max(y(3*end/5+1:4*end/5,:),-1),1) ; ...
                    min(max(y(4*end/5+1:end,:)-sigma*(R.*b),-mu),mu) ];
w = [0.7, 0.7, 0.5, 0.3];
A = @(x) [w(1)*Dx*x ; w(2)*Dy*x ; ...
          w(3)*DDz_new(x,numBands,numFrames) ; ...
          w(4)*Dt_new(x,numBands,numFrames); ...
          R.*STO(x)];
At = @(y) w(1)*Dxt*y(1:end/5,:) + w(2)*Dyt*y(end/5+1:2*end/5,:) + ...
          w(3)*DDzt_new(y(2*end/5+1:3*end/5,:),numBands,numFrames) + ...
          w(4)*Dtt_new(y(3*end/5+1:4*end/5,:),numBands,numFrames) + ...
          STO(R.*y(4*end/5+1:end,:));


%%  Determine timestep paramters 
opts = [];
%opts.L = .95/12;
opts.tol = .02;
%opts.tau = 0.5;  % 1.0 Works well for this problem when b has maximum size 10000
%opts.sigma = opts.L/opts.tau;
opts.maxIters = 2000;
opts.verbose = true;

[x,outs]= pdhg_adaptive(x0,y0,A,At,xProx, yProx,opts);

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







