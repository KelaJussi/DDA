function [path,full_loc,savename] = Make_DDA_save_folder_and_name(Input,save_type)
% Centralized function to make the right paths and filenames for saving DDA simulations.


% [sys] = Get_cpu_list_and_folder_hierarchy(true);
% Host-specific folders (optional, git-ignored: machine names and local paths)
if exist('Get_cpu_list_and_folder_hierarchy_TS', 'file') == 2
    [sys] = Get_cpu_list_and_folder_hierarchy_TS(false);
else
    sys = struct();
end
% Default data folder: DDA/Data next to DDA_Master.m (git-ignored)
data_folder = [fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Data') filesep];

% form string from array size parameters:
array_str = char(join(string(Input.arrayN),'_'));

% form a string variable from periodicities:
period_str = ['P' char(strrep(join(string(round(Input.p*1e9,1)),'_'),'.','p'))];

switch save_type
    %% Full results:
    case 'DDA_results'
        path = data_folder;
        
        if ~isfolder(path), mkdir(path); end
        
        savename=[Input.filename '.mat'];        
        full_loc=[path savename];
 
        %% NLO data:
    case {'SHG','SFG'}
        path=[sys.(Input.project).main_repo_folder sys.(Input.project).save_data_folder Input.m_name  '\'];
        mkdir(path)         
        
        file.name=[ save_type '_' Input.dimensionality '_'...
            Input.array_type '_'...
            'n_' strrep(num2str(Input.sample.n),'.','_') '_' ...
            period_str ...
            '_theta' strrep(num2str(Input.theta),'.','_') ...
            '_phi' strrep(num2str(Input.phi),'.','_') ...
            '_psi' strrep(num2str(Input.psi),'.','_') ...
            '_N'  array_str  ...
            '_tau'  strrep(num2str(round(Input.particle.tau(2)*1e15,2)),'.','_')     ...
            '_lspr' strrep(num2str(round(Input.particle.lambda_peak_vec(2)*1e9,1)),'.','_') ...
            '_cs' strrep(num2str(round(Input.particle.scatter_strength*1e8,2)),'.','_') ...
            '_'   num2str(Input.NLO_2D_limits(1)*1e9) ...
            '_to_' num2str(Input.NLO_2D_limits(2)*1e9) 'nm.mat'];

        savename=[file.name];
        full_loc=[path savename];
     
    
        %% Save plots:
    case '1D_spectra'
        path=[sys.(Input.project).main_repo_folder sys.(Input.project).save_fig_folder Input.m_name '\'];
        mkdir(path)
        
        savename=[Input.filename '.mat'];
        full_loc=[path savename];
        
    case '1D_plot_tmp'
        path=[Input.sys.save.tmp.data_folder{Input.sys.host_idn} Input.m_name '\'];
        mkdir(path)
        
        savename=[Input.filename '_' num2str(Input.II) '.mat'];
        full_loc=[path savename];
        
    case '1D_THG_spectra'
        path=[sys.(Input.project).main_repo_folder sys.(Input.project).save_data_folder Input.m_name  '\'];     
        mkdir(path)
        
        savename=[Input.filename '.mat'];
        full_loc=[path savename];
        
    case '2D_SFG_spectra'
        path=[sys.(Input.project).main_repo_folder  sys.(Input.project).save_data_folder Input.m_name '\'];
        mkdir(path)
        
        savename=['2DSFG_theta_' num2str(Input.theta) ...
                  '_phi_' num2str(Input.phi) ...
                  '_psi_' num2str(Input.psi) ...
                  '_chi_' Input.beta_struct '.mat'];    
        full_loc=[path savename];
    case '2D_SFG_spectra_fig'
        path=[sys.(Input.project).main_repo_folder  sys.(Input.project).save_data_folder Input.m_name '\'];
        mkdir(path)
        
        savename=['2DSFG_theta_' num2str(Input.theta) ...
                  '_phi_' num2str(Input.phi) ...
                  '_psi_' num2str(Input.psi) ...
                  '_chi_' Input.beta_struct '.fig'];    
        full_loc=[path savename];
        
    case '2D_DFG_spectra'
        path=[sys.(Input.project).main_repo_folder  sys.(Input.project).save_data_folder Input.m_name '\'];
        mkdir(path)
        
        savename=['2DDFG_theta_' num2str(Input.theta) ...
                  '_phi_' num2str(Input.phi) ...
                  '_psi_' num2str(Input.psi) ...
                  '_chi_' Input.beta_struct '.mat'];    
        full_loc=[path savename];
    case '2D_DFG_spectra_fig'
        path=[sys.(Input.project).main_repo_folder  sys.(Input.project).save_data_folder Input.m_name '\'];
        mkdir(path)
        
        savename=['2DSFG_theta_' num2str(Input.theta) ...
            '_phi_' num2str(Input.phi) ...
            '_psi_' num2str(Input.psi) ...
            '_chi_' Input.beta_struct '.fig'];
        full_loc=[path savename];
    otherwise
end

disp(['Save loc: ' full_loc])

end

