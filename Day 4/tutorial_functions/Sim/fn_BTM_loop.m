function [B11,B12,B21,B22] = fn_BTM_loop(u1,T1,p2,Qbar,Forcing,...
    Measurement,freqs)
    % Function that loops through several frequencies to compute the BTM.
    N = length(freqs);
    if N == 1
        B = Burner_Transfer_Matrix(u1,T1,p2,Qbar,Forcing,...
                Measurement, freqs,1);
            % Save values:
            B11 = B(1,1);
            B12 = B(1,2);
            B21 = B(2,1);
            B22 = B(2,2);
    else
        [B11,B12,B21,B22] = deal(zeros(N,1));
        parfor j = 1:N
            B = Burner_Transfer_Matrix(u1,T1,p2,Qbar,Forcing,...
                Measurement, freqs,j);
            % Save values:
            B11(j) = B(1,1);
            B12(j) = B(1,2);
            B21(j) = B(2,1);
            B22(j) = B(2,2);
        end
    end
end

function [B] = Burner_Transfer_Matrix(u1,T1,p2,Qbar,Forcing,Measurement,...
              freqs,j)
    % Function that computes the BTM at a single frequency
    Forcing.Frequency = freqs(j); % (Hz)
    Forcing.Speaker = 'ON';
    Forcing.Speaker_Downstream = 'OFF';
    % Upstream forcing ------------------------------------------------
    [Mean,~,Signals] = Run_Simulation(u1,T1,p2,Qbar,Forcing,Measurement);
    % Reshufle signals yo accomodate functions:
    S_up = Signals;
    S_dn = Signals;
    M_up = Measurement;
    M_dn = Measurement;
    % Keep only relevant Microphones:
    S_up.P_mic = S_up.P_mic(1:3,:);
    S_dn.P_mic = S_dn.P_mic(4:6,:);
    M_up.Mic_Pos = M_up.Mic_Pos(1:3);
    M_dn.Mic_Pos = M_dn.Mic_Pos(4:6);
    % MMM for Upstream forcing:
    % Multi-Microphone-Method:
    [Fu,Gu,~] = fn_MMM(Mean,S_up,M_up);
    [Fd,Gd,omega] = fn_MMM(Mean,S_dn,M_dn);
    % Signal Reconstruction:
    xu = 0;
    xd = 1e-5; % Position just after the jump
    [ph_u_uf,uh_u_uf] = fn_acoustic_field(Mean,Fu,Gu,omega,xu);
    [ph_d_uf,uh_d_uf] = fn_acoustic_field(Mean,Fd,Gd,omega,xd);
    % Downstream forcing ----------------------------------------------
    Forcing.Speaker = 'OFF';
    Forcing.Speaker_Downstream = 'ON';
    [Mean,~,Signals] = Run_Simulation(u1,T1,p2,Qbar,Forcing,Measurement);
    % Reshufle signals to accomodate functions:
    S_up = Signals;
    S_dn = Signals;
    M_up = Measurement;
    M_dn = Measurement;
    % Keep only relevant Microphones:
    S_up.P_mic = S_up.P_mic(1:3,:);
    S_dn.P_mic = S_dn.P_mic(4:6,:);
    M_up.Mic_Pos = M_up.Mic_Pos(1:3);
    M_dn.Mic_Pos = M_dn.Mic_Pos(4:6);
    % MMM for Upstream forcing:
    % Multi-Microphone-Method:
    [Fu,Gu,~] = fn_MMM(Mean,S_up,M_up);
    [Fd,Gd,omega] = fn_MMM(Mean,S_dn,M_dn);
    % Signal Reconstruction:
    [ph_u_df,uh_u_df] = fn_acoustic_field(Mean,Fu,Gu,omega,xu);
    [ph_d_df,uh_d_df] = fn_acoustic_field(Mean,Fd,Gd,omega,xd);
    % Compute matrix:
    B = fn_BTM(Mean, ph_d_uf, uh_d_uf, ph_d_df, uh_d_df,...
                     ph_u_uf, uh_u_uf, ph_u_df, uh_u_df);
end