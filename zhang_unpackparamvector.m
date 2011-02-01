% decode the optimisation parameter vector x and store all the relevant
% properties inside the calibration structure
function Calibration = zhang_unpackparamvector( calib, x )

Calibration = zhang_load( calib );

if Calibration.square_aspect
    focal = x([1,1]);
    skewIdx = 2;
else
    focal = x([1,2]);
    skewIdx = 3;
end
if Calibration.zero_skew
    skew = 0;
    principleIdx = skewIdx;
else
    skew = x(skewIdx);
    principleIdx = skewIdx + 1;
end
principle = x([principleIdx,principleIdx+1]);
radialIdx = principleIdx + 2;
if Calibration.first_order
    Calibration.k1 = x(radialIdx);
    rtIdx = radialIdx + 1;
else
    Calibration.k1 = x(radialIdx);
    Calibration.k2 = x(radialIdx+1);
    rtIdx = radialIdx + 2;
end

Calibration.A = [ focal(1)   skew       principle(1)
                  0          focal(2)   principle(2)
                  0          0          1             ];

active = [Calibration.Images.Active];
nActive = nnz( active );
nPlanes = (length(x)-(rtIdx-1)) / 6;
if nPlanes ~= nActive
    error( 'x vector length incompatible with active image set' );
end


nthActiveImage = 0;
for i = 1:length( Calibration.Images )
    if Calibration.Images(i).Active
        nthActiveImage = nthActiveImage + 1;
        first = (nthActiveImage-1)*6 + rtIdx;
        Calibration.Images(i).R = rodrigues_vec2mat( x(first:first+2) );
        Calibration.Images(i).T = x(first+3:first+5);
        % [brad] bugfix: Dec 2010: added missing set of Rv parameter
        Calibration.Images(i).Rv = rodrigues_mat2vec( Calibration.Images(i).R );
     else
        Calibration.Images(i).R = [];
        Calibration.Images(i).T = [];
        Calibration.Images(i).Rv = [];
    end
end


% output
if ischar( calib )
    disp( 'saving to disk' );
    save( calib, 'Calibration', '-append' );
    disp( 'done' );
end

    
