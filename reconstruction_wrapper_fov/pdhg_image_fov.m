function [x,outs ] = pdhg_image_fov( R, b, mu, order, res_prev, topLeft, botRight, record_ds, fovParams, guess )

rows = order.n; 
cols = order.n;
order_prev = createOrderingData(res_prev,'full');

%%  Allocate memory
rows = order.n1;
cols = order.n2;
x0=zeros(cols,1);
y0=zeros(2*cols,1);

%% Get Derivative Matrices
[Dx,Dy] = createDifferenceOperators_fov(order,order_prev,topLeft,botRight);
Dxt = Dx';
Dyt = Dy';


b = b*100/max(max(abs(b)));
rhs = mu*fSTO_RLE(R.*b,fovParams);

xProx =@(x,tau) fSTO_RLE( ...
       (1/tau+mu*R).^-1.*STO(UndoRLE(rhs+x/tau,record_ds,fovParams)),fovParams);
yProx = @(y,sigma) projectInf(y(1:end/2,:),y(end/2+1:end,:) );
A = @(x) [Dx*x ; Dy*x];
At = @(y) Dxt*y(1:cols,:)+Dyt*y(cols+1:2*cols,:);

%%  Determine timestep paramters 
opts = [];
% opts.L = 0.95/8;
opts.tol = .01;
% opts.tau = 2;  % 1.0 Works well for this problem when b has maximum size 10000
% opts.sigma = opts.L/opts.tau;
opts.maxIters = 500;


if exist('guess','var')
x0 = guess;
y0 = [sign(Dx*x0) ; sign(Dy*x0)];
end

[x,outs]= pdhg_adaptive(x0,y0,A,At,xProx, yProx,opts);

  


% %en = 0;
% rp = 1;
% rd = 1;
% 
% 
% rhs = mu*STO(R.*b);
% 
% iter = 1;
% while max(rp(end)/rp(1),rd(end)/rd(1))>0.01 && iter<500
% %while iter<50 
% 
%     u0 = u;
%     x0 = x;
%     y0 = y;
%     
%     uhat = u-tau*(Dxt*x+Dyt*y);
%     
%     u = rhs+uhat/tau;
%     u = STO(u);
%     u = u./(mu*R+1/tau);
%     u = STO(u);
%     
%   
%     uh = u+(u-u0);
%   
%     
%     x = x+sigma*Dx*uh;
%     y = y+sigma*Dy*uh;
%     normalize = sqrt(x.*x+y.*y);
%     normalize = max(normalize,1);
%     x = x./normalize;
%     y = y./normalize;
%     
%    
%     rp(iter) = norm(u-u0,'fro')/tau;
%     rd(iter) = sqrt(norm(x-x0,'fro')^2+norm(y-y0,'fro')^2)/sigma/2.0;
%     %en(iter) = sum(sum( sqrt((Dx*u).^2+(Dy*u).^2)))+mu/2*norm(R.*hadamardTransform(B.*u)/scale-s,'fro')^2;
%    
%    
%     disp([num2str(iter) ' : ' num2str(rp(end))]);
%     
%     iter = iter+1;
%     
%     
% end
% 
% image = nestedVectorToImage(u, data);
%     
% return
