function [loop_counter]=Display_calc_times_and_convergence_specs(Process_tag,loop_counter,I_wl,modality_timer_tag,flag,Res,t0)
% Displays DDA calculation times and informs if problems in convergence.

if nargin < 5 || isnan(Res)
    switch Process_tag
        case {'sfg','SFG'}
            %% 5a. display calculation times of the first 10 loops:
            if (mod(loop_counter,10) == 0 && loop_counter < 12)
                disp([num2str(loop_counter) ' / ' num2str(round(length(I_wl)*length(I_wl)/2)) modality_timer_tag ...
                    num2str(toc,4) ' s.']);
            end
            
            % 5b. display calculation times of every 1000th loop:
            if mod(loop_counter,250) == 0
                disp([num2str(loop_counter) ' / ' num2str(round(length(I_wl)*length(I_wl)/2)) modality_timer_tag ...
                    num2str(toc,4) ' s.']);
            end
            loop_counter=loop_counter+1;
        
        case {'dfg','DFG'}
            %% 5a. display calculation times of the first 10 loops:
            if (mod(loop_counter,10) == 0 && loop_counter < 12)
                disp([num2str(loop_counter) ' / ' num2str(round(length(I_wl)*length(I_wl)/2)) modality_timer_tag ...
                    num2str(toc,4) ' s.']);
            end
            
            % 5b. display calculation times of every 1000th loop:
            if mod(loop_counter,1000) == 0
                disp([num2str(loop_counter) ' / ' num2str(round(length(I_wl)*length(I_wl)/2)) modality_timer_tag ...
                    num2str(toc,4) ' s.']);
            end
            loop_counter=loop_counter+1;
            
        case 'SHG'
            if mod(loop_counter,10) == 0
                disp([num2str(loop_counter) ' / ' num2str(length(I_wl)) modality_timer_tag ...
                    num2str(toc,4) ' s.']);
            end
            loop_counter=loop_counter+1;
            
        case 'fwm'
            if mod(loop_counter,10) == 0
                disp([num2str(loop_counter) ' / ' num2str(length(I_wl)) modality_timer_tag ...
                    num2str(toc(t0),4) ' s.']);
            end
            loop_counter=loop_counter+1;
        case 'fwm100'
            if mod(loop_counter,100) == 0
                disp([num2str(loop_counter) ' / ' num2str(length(I_wl)) modality_timer_tag ...
                    num2str(toc(t0),4) ' s.']);
            end
            loop_counter=loop_counter+1;
    end
else
    if flag ~= 0
        Convergence_counter = flag;
        disp([' Error in convergence: Flag = ' num2str(Convergence_counter) '. Residual = ' num2str(Res)])
    end
end




end