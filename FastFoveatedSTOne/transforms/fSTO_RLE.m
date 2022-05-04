%%  
%  This function performs the STO transform by calling the c implementation
%  of STO_fast.  This method can be called on either a vector or a matrix.
%  When called on a matrix, it performs a separate transform on each
%  column.


function [ out ] = fSTO_RLE( in, fovParams)

assert(ndims(in)<3, 'Input must be a vector or matrix');
numRLE_px = fovParams(1)^2 - length(fovParams(3:end)) + ...
            length(fovParams(3:end))/(fovParams(1)/fovParams(2))^2;

%  Handle case that 'in' is a row vector
if min(size(in))==1
    out = in(1:numRLE_px);
    %out = in;
    % This line forces matlab NOT to use a lazy copy.  This is important
    % because otherwise a call to "STO_fast" will overwrite the input
    out(1)=out(1);
    fSTO_fast_RLE(in,fovParams);
    out = in(1:numRLE_px);

else
    for c=1:size(in,2)
        slice = in(:,c);
        fSTO_fast_RLE(slice,fovParams);
        out(:,c) = slice(1:numRLE_px,:);
    end
end
