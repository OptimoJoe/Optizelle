% Check that we have an equality constrained state 
function checkEqualityConstrainedStateT(name,value)
    % Set the error message
    err = sprintf( ...
        'The %s argument must have type EqualityConstrained.State.t.', ...
        name);

    % Check the unconstrained values
    try
        checkUnconstrainedStateT(name,value);
    catch
        error(err);
    end

    % Check for the appropriate fields 
    if ~(checkFields({ ...
        'y', ...
        'dy', ...
        'zeta', ...
        'eta0', ...
        'rho', ...
        'rho_old', ...
        'rho_bar', ...
        'eps_constr', ...
        'xi_qn', ...
        'xi_pg', ...
        'xi_proj', ...
        'xi_tang', ...
        'xi_lmh', ...
        'xi_lmg', ...
        'xi_4', ...
        'rpred', ...
        'PSchur_left_type', ...
        'PSchur_right_type', ...
        'augsys_iter_max', ...
        'augsys_rst_freq', ...
        'g_x', ...
        'norm_gxtyp', ...
        'gpxdxn_p_gx', ...
        'gpxdxt', ...
        'norm_gpxdxnpgx', ...
        'dx_n', ...
        'dx_ncp', ...
        'dx_t', ...
        'dx_t_uncorrected', ...
        'dx_tcp_uncorrected', ...
        'H_dxn', ...
        'W_gradpHdxn', ...
        'H_dxtuncorrected', ...
        'g_diag'}, ...
        value))
        error(err);
    end
end