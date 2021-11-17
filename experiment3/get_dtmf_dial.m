function dial = get_dtmf_dial(fr_calc, fc_calc)
%GET DIEAL Summary of this function goes here
%   Detailed explanation goes here

    % Estimate frequency of the approximated DTFT of the DTMF signal
    dtmf_keys = ['1','2','3','A';
                 '4','5','6','B';
                 '7','8','9','C';
                 '*','0','#','D'];
    f_row = [697 770 852 941];
    f_col = [1209 1336 1477 1633];
    
    [~, pos_row] = min(abs(f_row - fr_calc)); % Get the position of the nearest row frequency 
    [~, pos_col] = min(abs(f_col- fc_calc)); % Get the position of the nearest collumn frequency
    
    dial = dtmf_keys(pos_row, pos_col);  % extract the dial of the input signal
end

