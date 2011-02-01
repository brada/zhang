% after running zhang_calibrate_optimise, we may want to impose some
% nonlinear equality constraints. This function will try to do that. Note
% however that this will usually take MUCH longer to compute, and the
% reprojection errors will get a lot worse (at first) while it tries to
% satisfy the constraint that was probably not satisfied by
% calibrate_optimise
function [Calibration,x] = zhang_calibrate_optimise_nonlcon( calib, debug )

Calibration = zhang_load( calib );
x0 = Calibration.x;

% check for constraint satisfaction before running fmincon
[~,ceq] = zhang_constraint_eq( Calibration, x0 );
if norm(ceq) == 0
    disp( 'Constraints are already satisfied exactly' );
elseif norm(ceq) < 1e-6
    disp( 'Constraints are already satisfied closely' );
else
    
    if debug
        plotfuncs   = {@optimplotx,@optimplotfirstorderopt,@optimplotresnorm,@optimplotstepsize};
        display     = 'iter';
        diagnostics = 'on';
        tolX        = 1e-6;
        tolFun      = 1e-6;
    else
        diagnostics = 'off';
        plotfuncs   = [];
        display     = 'final';
        tolX        = 1e-10;
        tolFun      = 1e-10;
    end
    
    % default trust-region-reflective algorithm does not work well
    % (possibly because I don't supply a gradient of obj_fun_scalar).
    % interior-point works well even if the constraints are not very
    % closely satisfied by the initial guess
    options = optimset( ...
        'Algorithm', 'interior-point', ...
        'Diagnostics', diagnostics, ...
        'PlotFcns', plotfuncs, ...
        'TolX', tolX, ...
        'TolFun', tolFun, ...
        'MaxIter', 100, ...
        'MaxFunEvals', 1000000*length(Calibration.x0), ...
        'TypicalX', Calibration.x, ...
        'SubproblemAlgorithm', 'cg', ...
        'Display', display );
                
    disp( 'Running constrained optimisation...' );
    obj_fcn = @(x) zhang_reprojectionerror( Calibration, x );
    % this matches the f(x) reported by lsqnonlin in
    % zhang_calibrate_optimise. Note that mean(obj_fcn(x)) and
    % norm(obj_fcn(x)) do NOT match
    obj_fcn_scalar = @(x) sum( obj_fcn(x).^2 );
    nonlcon = @(x) zhang_constraint_eq( Calibration, x );
    [x,~,~,output,~,grad] = fmincon( obj_fcn_scalar, Calibration.x, [],[],[],[],[],[], nonlcon, options );
    
    if debug
        disp( 'Constraint violation of unconstrained optimisation:' );
        disp( ceq );
        disp( 'Reprojection errors of unconstrained optimisation:' );
        disp( Calibration.mean_reprojection_error_pixels );
    end
       
    
    % update calibration with new parameters
    Calibration.x = x;
    Calibration.nonlconGradient = grad;
    Calibration.nonlconOutput = output;
    if debug
        disp( 'grad:' );
        disp( grad' );
        disp( 'output:' );
        disp( output );
    end

    % note that fmincon takes a scalar function, not a vector like
    % lsqnonlin, so we are minimising just the mean of parmeter vector x,
    % not x itself. This also means that from fmincon we get an output
    % gradient, not a jacobian. Lack of the jacobian means we cannot now
    % compute new confidence intervals, so Calibration.confidence_interval
    % will reflect the confidence of the solution after
    % zhang_calibrate_optimise, NOT after zhang_calibrate_optimise_nonlcon
    
    % apply the current solution to all the image extrinsic parameters and
    % project the points according to this solution
    Calibration = zhang_unpackparamvector( Calibration, x );
    Calibration = zhang_projectpoints( Calibration );

    % get the mean reprojection error across all images
    mean_re_ndc = [Calibration.Images(:).mean_reprojection_error_NDC];
    Calibration.mean_reprojection_error_NDC = mean_re_ndc;
    Calibration.mean_reprojection_error_pixels = mean_re_ndc / Calibration.N(1,1);
    
    if debug
        disp( 'Reprojection errors of constrained optimisation:' ); 
        disp( Calibration.mean_reprojection_error_pixels );
    end
    
end
