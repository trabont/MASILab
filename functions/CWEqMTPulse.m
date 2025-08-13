function w1e = CWEqMTPulse(w1,t,teq)
% this function calculates the continuous wave
% equilivant (in terms of RF power) of saturation MT pulse
% w1 : pulse amplitude (rad/s)
% t  : corresponding time points of pulse (s)
% teq: length of equilivant pulse (s)
% w1e : continuous wave equilivant w1 (rad/s)

% Loop for each RF pulse supplied (each column of w1 and t)
for ii = 1:size(w1,2)
    t_in = 6.4e-6/2:6.4e-6:t(end,ii);
    B1MTi = interp1(t(:,ii),w1(:,ii),t_in,'pchip');
    % Integrate w1^2
    w12int = trapz(t_in,B1MTi.^2);
    
    % Determine cw equilivant square pulse 
    w1e(ii) = sqrt(w12int/teq(ii));
end
