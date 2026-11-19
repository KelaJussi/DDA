function [] = DDA_simulation_error(error_number)
% List of obvious errors.

if error_number == 1
    disp('Current computer unknown.')
    disp('Please update the details into')
    disp('Get_cpu_list_and_folder_hierarchy.m')
    error('myApp:argChk','FIX: Set path -> add with sub-folders to: /General_codes/')
end

end

