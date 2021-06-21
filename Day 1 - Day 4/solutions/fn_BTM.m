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
    y = [ph_d_uf/(rho2*c2);uh_d_uf;ph_d_df/(rho2*c2);uh_d_df];
    M = [ph_u_uf/(rho1*c1) uh_u_uf 0 0; ...
         0 0 ph_u_uf/(rho1*c1) uh_u_uf; ...
         ph_u_df/(rho1*c1) uh_u_df 0 0; ...
         0 0 ph_u_df/(rho1*c1) uh_u_df]; 
    b = M\y;
    
    B = [b(1) b(2); b(3) b(4)];
    % Sanity check, Remember that B must be a 2 x 2 matrix!
    sz = size(B);
    if sz(1) ~= 2 && sz(2) ~= 2
        error('Matrix B is not 2 x 2')
    end
end