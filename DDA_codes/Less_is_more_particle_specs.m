function [particle] = Less_is_more_particle_specs(particle_type)
% Use this centralized function to ensure that similar parameters are used
% in all simulations.

particle.particle_type=particle_type;

% scatter_strength_A1R=0.18e-5/(4*pi);
scatter_strength_A1R=0.20e-5/(4*pi);
% lambda_peak_vec_A1R=1e-9*[1550 1094 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
% lambda_peak_vec_A1R=1e-9*[1550 1165 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
lambda_peak_vec_A1R=1e-9*[1550 1145 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
lambda_peak_vec_A1R=1e-9*[1550 1100 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
% lambda_peak_vec_A1R=1e-9*[1550 1060 800]; % plasmon resonances lambda_x, lambda_y, lambda_z

% lambda_peak_vec_A1R=1e-9*[1550 1145 800]; % plasmon resonances lambda_x, lambda_y, lambda_z

% % lambda_peak_vec_A1R=1e-9*[1550 1135 800]; % plasmon resonances lambda_x, lambda_y, lambda_z

tau_A1R=[2.5 1.9 1.85]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
    
if nnz(strfind(particle_type,'L175'))
    %% These single particle values are gotten by comparing single particle results:
    % L = 175 nm particles)
    % Input.lambda_peak_vec=1e-9*[832 832 832]; % plasmon resonances lambda_x, lambda_y, lambda_z
    % Input.particle.tau=[0.9750  0.9750  0.9750]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
    
    
    % L = 175 nm particles), lspr peak from measurements ad hoc)
%     particle.lambda_peak_vec=1e-9*[1610 1050 1050]; % plasmon resonances lambda_x, lambda_y, lambda_z
    %% 21.8. efforts to have the measured and simulated coordinate systems coincide...
    % meas: x-pol gives shorter wl LSPR , y-pol longer
    particle.lambda_peak_vec=1e-9*[1052 600 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
    particle.tau=[1.9 1.92 1.85]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
     particle.scatter_strength=0.15e-5/(4*pi); % A0=0.3e-5 1Dz: 0.05e-5!!! THIS IS HUGELY IMPORTANT !!

     particle.scatter_strength_range=[0.05 1];
     
elseif nnz(strfind(particle_type,'1a'))
    % L = 275 nm particles)
    %     particle.lambda_peak_vec=lambda_peak_vec_A1R; % plasmon resonances lambda_x, lambda_y, lambda_z
    %     particle.tau=tau_A1R;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
    %     particle.scatter_strength=scatter_strength_A1R; % A0=0.3e-5 1Dz: 0.05e-5!!! THIS IS HUGELY IMPORTANT !!
    particle.tau=[2.5 1.6 1.85]*1e-15;
    particle.scatter_strength=0.125e-5/(4*pi); % 0.13 was pretty close for N=100
%     particle.scatter_strength=0.20e-5/(4*pi);
    particle.lambda_peak_vec=1e-9*[1550 1095 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
    particle.scatter_strength_range=[0.05 1];
    
    %% extinction data filename:
    particle.ext_data_filename.y='lowerR1a_90deg.dat';
    particle.ext_data_filename.x='lowerR1a_0deg.dat';
%     particle.ext_data_filename.y='lowerL1a_90deg.dat';
%     particle.ext_data_filename.x='lowerL1a_0deg.dat';
    
elseif nnz(strfind(particle_type,'1b'))
    % L = 275 nm particles)
%     particle.lambda_peak_vec=lambda_peak_vec_A1R; % plasmon resonances lambda_x, lambda_y, lambda_z
%     particle.tau=tau_A1R;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
%     particle.scatter_strength=scatter_strength_A1R; % A0=0.3e-5 1Dz: 0.05e-5!!! THIS IS HUGELY IMPORTANT !!

    particle.tau=[2.5 1.75 1.85]*1e-15;
    particle.scatter_strength=0.13e-5/(4*pi); % 0.2
    particle.lambda_peak_vec=1e-9*[1550 1137 800]; % plasmon resonances lambda_x, lambda_y, lambda_z

%     particle.lambda_peak_vec=1e-9*[1550 1146 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
%     particle.tau=[2.5 1.95 1.85]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
    
%     particle.scatter_strength=scatter_strength_A1R; % A0=0.3e-5 1Dz: 0.05e-5!!! THIS IS HUGELY IMPORTANT !!    
    
    %% extinction data filename:
    particle.ext_data_filename.y='lowerR1b_90deg.dat';
    particle.ext_data_filename.x='lowerR1b_0deg.dat';
    
elseif nnz(strfind(particle_type,'1c'))
    % L = 275 nm particles)
    particle.tau=[2.5 0.9 1.85]*1e-15;
    particle.scatter_strength=0.1e-5/(4*pi); % 0.13e-5/(4*pi) 0.08
	%     particle.lambda_peak_vec=1e-9*[1550 1160 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
    particle.lambda_peak_vec=1e-9*[1550 1140 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
    %% extinction data filename:
    particle.ext_data_filename.y='lowerR1c_90deg.dat';
    particle.ext_data_filename.x='lowerR1c_0deg.dat';

elseif nnz(strfind(particle_type,'1d'))
    % L = 275 nm particles)
    particle.tau=[2.5 1.7 1.85]*1e-15;
    particle.scatter_strength=0.18e-5/(4*pi); % 0.2
    particle.lambda_peak_vec=1e-9*[1550 1143 800]; % plasmon resonances lambda_x, lambda_y, lambda_z
    %         particle.lambda_peak_vec=1e-9*[1550 1137 800]; % plasmon resonances lambda_x, lambda_y, lambda_z

    %% extinction data filename:
    particle.ext_data_filename.y='lowerR1d_90deg.dat';
    particle.ext_data_filename.x='lowerR1d_0deg.dat';
    
elseif nnz(strfind(particle_type,'L250'))
        % L = 250 nm particles)
%     Input.lambda_peak_vec=1e-9*[1230 1230 1230]; % plasmon resonances lambda_x, lambda_y, lambda_z
%     Input.particle.tau=[1.6195 1.6195 1.6195]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15

elseif nnz(strfind(particle_type,'C1v')) 
%     save beta_tmp
    
    %% 21.8. trial
    C=1;
    particle.beta_omega3(2,2,2)=0.1*C;                       % yyy
    particle.beta_omega3(2,1,1)=C;                       % -yxx
    particle.beta_omega3(1,1,2)=0.001*C;                       % -xxy
    particle.beta_omega3(1,2,1)=beta_omega3(1,1,2);      % -xyx
    
    
else
        % L = 175 nm particles), lspr peak from measurements ad hoc)
    particle.lambda_peak_vec=1e-9*[1610 1060 1060]; % plasmon resonances lambda_x, lambda_y, lambda_z        
    particle.tau=[1.15 1.15 1.15]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
     
    particle.scatter_strength=0.26e-5/(4*pi); % A0=0.3e-5 1Dz: 0.05e-5!!! THIS IS HUGELY IMPORTANT !!
  
end


end

