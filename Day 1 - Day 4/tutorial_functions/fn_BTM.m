function [B] =  fn_BTM(Mean, ph_d_uf, uh_d_uf, ph_d_df, uh_d_df,...
                             ph_u_uf, uh_u_uf, ph_u_df, uh_u_df)
    % Function that computes the Burner Transfer Matrix.
    
    % Read the mean flow properties:
    rho1 = Mean.rho1;
    rho2 = Mean.rho2;
    c1 = Mean.c1;
    c2 = Mean.c2;
    % ====================================================================
    % Using Eq. (4) compute the matrix B:



    % Sanity check, Remember that B must be a 2 x 2 matrix!
    sz = size(B);
    if sz(1) ~= 2 && sz(2) ~= 2
        error('Matrix B is not 2 x 2')
    end
end
