%%  
%  This function undoes the RLE by calling the c implementation of
%  STO_fast.  This method can be called on either a vector or a matrix.
%  When called on a matrix, it acts separately on each column.

function [ out ] = UndoRLE(x,record_ds,fovParams)

assert(ndims(x)<3, 'Input must be a vector or matrix');
jump_sqrt = fovParams(1)/fovParams(2);
jump = jump_sqrt*jump_sqrt;
numOutPx = uint32(sum(record_ds==1)*jump+sum(record_ds==0));

%  Handle case that 'in' is a row vector
record_ds = logical(record_ds);
jump_sqrt = uint32(jump_sqrt);
if min(size(x))==1
    out = UndoRLE_fast(x,record_ds,jump_sqrt,numOutPx);
else
    out = zeros(numOutPx,size(x,2));
    for i = 1:size(x,2)
        out(:,i) = UndoRLE_fast(x(:,i),record_ds,jump_sqrt,numOutPx);
    end
end
