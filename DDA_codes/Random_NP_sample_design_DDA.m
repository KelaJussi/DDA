function [Xout,Yout]=Random_NP_sample_design_DDA(Input)
% Function to make an array consisting of NxN randomly positioned
% nanoparticles.                                 - mjh 26.10.2017
% 
% Idea is that, when the distances between particles are randomized, no
% collective responses occur in the array, and the overall response is
% close to the response of an individual particle. This way, measuring e.g.
% extinction spectra etc. is a lot easier.   
%
% [ Xout, Yout ] = x- and y-coordinates for the particles in randomized
%                  array of given size 

%% Define array configuration
Array.N=Input.arrayN; %100;                        % N^2 = total number of particles in the array 
Array.size=Input.array_size; %200e-6;                  % physical array size
Array.min_particle_distance=800e-9; % minimum particle separation

plot_handle=false;                  % handle to control resulting array scatter plot

% first trial for the randomized position vector
scaler=1.5;                         % extra factor to ensure that we can still have Array.N^2 particles in array after deletion 
Array.X=(Array.size)*rand(Array.N^2*scaler,1)-Array.size/2;   
Array.Y=(Array.size)*rand(Array.N^2*scaler,1)-Array.size/2;

%% calculate the distances between checked particles:
Xtmp=Array.X;
Ytmp=Array.Y;

for jj=1:length(Array.X) % ii:(ii+x_sweep_stride)
    if length(Xtmp)<jj
        % calculate the distances between checked particles:
        Dist_tmp=sqrt((Xtmp-Xtmp(end)).^2+(Ytmp-Ytmp(end)).^2);
        Dist_tmp=Dist_tmp(Dist_tmp>0);
        I=find(Dist_tmp < Array.min_particle_distance & Dist_tmp>0);
        
        % delete particles which are too close to other particles:
        if nnz(I)
            Xtmp(I)=[];
            Ytmp(I)=[];
        end
        break
    end
    
    % calculate the distances between checked particles:
    Dist_tmp=sqrt((Xtmp-Xtmp(jj)).^2+(Ytmp-Ytmp(jj)).^2);
    I=find(Dist_tmp < Array.min_particle_distance & Dist_tmp>0);
    
    % delete particles which are too close to other particles:
    if nnz(I)
        Xtmp(I)=[];
        Ytmp(I)=[];
    end
end

%% take the first N^2 particles, if enough exist. 
if length(Xtmp) < Array.N.^2
    disp(['Problem in fitting enough particles into the array.'])
    disp(['Try again by decreasing minimum particle separation or scaler..'])
    
    Xout=[];
    Yout=[];
    return
else
    Xout=Xtmp(1:Array.N.^2);
    Yout=Ytmp(1:Array.N.^2);
end

disp(['Trial particle N^2 = ' num2str(length(Array.X)) ... 
      ', randomized array N^2 = ' num2str(length(Xtmp)) ]) 
      

if plot_handle
    %% Plotting:
    sub_scaling=[0.08 0.08];
    
    figure
    subplot_tight(1,2,1,sub_scaling)
    hold on
    scatter(Xout,Yout,'x')
    for aa=1:length(Xtmp)
        circle(Xout(aa),Yout(aa),Array.min_particle_distance)
    end
    axis equal
    
    subplot_tight(1,2,2,sub_scaling)
    scatter(Xout,Yout,'x')
    axis equal
    
    set(gcf,'position',[132         147        1331         768])
end

end