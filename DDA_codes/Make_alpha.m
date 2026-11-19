function [alpha_diag]=Make_alpha(Array,Input,kk)

alpha_diag=zeros(Input.N,1,'like',Array.r_i);
switch Input.dimensionality
    %% make N x1 vector too out of the diagonal components:
    case '1Dx'
        alpha_diag(1:end)=Array.alpha(:,1,1,kk);
    case '1Dy'
        alpha_diag(1:end)=Array.alpha(:,2,2,kk);
    case '1Dz'
        alpha_diag(1:end)=Array.alpha(:,3,3,kk);
    case '3D'
        alpha_diag=zeros(3*Input.N,1,'like',Array.r_i);
        %% make 3Nx1 vector out of the diagonal components:
        for ii=1:3
            alpha_diag(ii:3:end)=Array.alpha(:,ii,ii,kk);
        end
    case '2D'
        alpha_diag=zeros(2*Input.N,1,'like',Array.r_i);
        for ii=1:2
            alpha_diag(ii:2:end)=Array.alpha(:,ii,ii,kk);
        end
    case 'cyl_TM'
        alpha_diag(1:end)=Array.alpha(:,kk);
    otherwise
        alpha_diag(1:end)=Array.alpha(:,1,1,kk);
end

% ad hoc (but working) way of having the FCC array for Czaplicki's Less Is More -stuff:
switch Input.array_type
    case 'fcc_square'
        for xx=1:3
            alpha_diag(xx:3:end)=Array.checkerboard(1:end).'.*alpha_diag(xx:3:end);
        end
    otherwise
end

% ad hoc (but working) way of having the hexagonal array:
switch Input.array_type
    case 'hexagonal'
        for xx=1:3
            alpha_diag(xx:3:end)=Array.checkerboard(1:end).'.*alpha_diag(xx:3:end);
        end
    otherwise
end

end