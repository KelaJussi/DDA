function [Array,Input]=Make_generic_Array(Input)
% makes the Array struct, forgets z-components in positions & vectors 
% outputs: r_i = position vector Nx2
 
%% switch to control which sample type is used:
switch Input.sample.type
    case 'lorentzian Au'          
        if 0
            [Array,Input]=Make_3D_Array_sub_JTK(Input);
        else
            [Array,Input]=Make_3D_Array_sub(Input);  
        end
         
    case 'oligomer'
        [Array,Input]=Make_3D_Oligomer_sub(Input);   
    case 'stiched_array'
%         save tmp_stitch
         [Array,Input]=Make_3D_Array_sub(Input);          
    case 'random_surface'
%         [Array,Input]=Make_3D_Array_sub(Input);  
         [Array,Input]=Make_Random_Array_sub(Input);   
    otherwise
        [Array,Input]=Make_3D_Array_sub(Input);          
end

Input.N=size(Array.r_i,1); % number of particles in array
%disp(['Array formed in ' num2str(toc-Input.toc1) ' s'])

end

function [Array,Input]=Make_3D_Array_sub(Input)   

%% nested function to find the correct hex array size using fzero
    function x=find_root_for_correct_hex_array(N)
        x=(Input.arrayN^2-2*N.^2);
    end

switch  Input.array_type

    case 'square'      
        
        if length(Input.arrayN) >2
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
            array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
            array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);            
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
            r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
            r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
            
            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;       
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end              
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
                 
        % build Nx2 ri as output position vector
        for ii=1:3 % 3
             Array.r_i(:,ii)=single(R_i{ii}(1:end));
        end
        
        switch Input.sample.type
            case 'stiched_array'
                % shift the quadrants of the array to mimick stitching
                % error:
                
                Ix=find(Array.r_i(:,1)>array_max_xlength/2);    % x sticthing:
                Iy=find(Array.r_i(:,2)>array_max_ylength/2);    % y sticthing:
                Array.r_i(Ix,1)=Input.sample.sticthing(1)+Array.r_i(Ix,1);
                Array.r_i(Ix,2)=Input.sample.sticthing(2)+Array.r_i(Ix,2);             
                Array.r_i(Iy,2)=Input.sample.sticthing(2)+Array.r_i(Iy,2);
                
%                 %% plot the array for testing:
%                 figure
%                 scatter(Array.r_i(:,1),Array.r_i(:,2),'x')
%                 axis square
%                 set(gcf,'position',[972  569  894  470])
            
            otherwise
        end
        
        
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            max_size(ii)=Input.p(ii)*(Input.Nmax(ii)-1);
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
        end       

        %         %% plot the array for testing:
        %         figure
        %         scatter(Array.r_i_conv(:,1),Array.r_i_conv(:,2),'x')
        %         set(gcf,'position',[139   156  1594   923])
        
    case 'shifted_diagonals'
        if length(Input.arrayN) >2
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
            array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
            array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);            
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
            r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
            r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
            
            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;       
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end   



        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
         x = mod([1:length(r_xi)],2)*2-1;
         y = mod([1:length(r_yi)],2)*2-1;
        checkerboard = (y'*x) > 0;
        R_i{1}(checkerboard) = R_i{1}(checkerboard)-Input.diag_shift;
        R_i{2}(checkerboard) = R_i{2}(checkerboard)-Input.diag_shift;
        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end

               
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            %max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end       

                %% plot the array for testing:
                figure
                scatter(Array.r_i(:,1),Array.r_i(:,2),'x')
                set(gcf,'position',[139   156  1594   923])
    
    case '1Dsquare'      
        
        if length(Input.arrayN) >2
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
            array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
            array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);            
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
            r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
            r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
            
            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;       
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end              
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
              
        %% add slight randomness to the the other (x) direction coordinates here:     
        if 0
        for ii=1:Input.Nmax(1)
            % 8.5. it seems that the randomness wasn't enough
             R_i{1}(ii,:)=R_i{1}(ii,:)+0.5*Input.p(1)*rand([1 1]);
%              R_i{1}(ii,:)=R_i{1}(ii,:)+*Input.p(1)*rand([1 1]);
        end
        end
%         R_i{1}=R_i{1}+Input.p(1)*rand(size(R_i{1}));
            
%         save tmp_1D
        % build Nx2 ri as output position vector
        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        
%         figure; scatter(Array.r_i(:,1),Array.r_i(:,2)); axis square
        
        switch Input.sample.type
            case 'stiched_array'
                % shift the quadrants of the array to mimick stitching
                % error:                
                Ix=find(Array.r_i(:,1)>array_max_xlength/2);    % x sticthing:
                Iy=find(Array.r_i(:,2)>array_max_ylength/2);    % y sticthing:
                Array.r_i(Ix,1)=Input.sample.sticthing(1)+Array.r_i(Ix,1);
                Array.r_i(Ix,2)=Input.sample.sticthing(2)+Array.r_i(Ix,2);             
                Array.r_i(Iy,2)=Input.sample.sticthing(2)+Array.r_i(Iy,2);
        
            otherwise
        end
               
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            %max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end       

        %         %% plot the array for testing:
        %         figure
        %         scatter(Array.r_i_conv(:,1),Array.r_i_conv(:,2),'x')
        %         set(gcf,'position',[139   156  1594   923])
    
    case '1Dsquare_halfshifts'      
        
        if length(Input.arrayN) >2
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
            array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
            array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);            
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
            r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
            r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
            
            for ll = 1:(Input.array_shifts_y)
                
                 r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)=r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)+Input.p(2)/2;
            end
            for ll = 1:(Input.array_shifts_x)
                
                 r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)=r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)+Input.p(2)/2;
            end

            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
            
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;       
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end              
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
              
        %% add slight randomness to the the other (x) direction coordinates here:     
        if 0
        for ii=1:Input.Nmax(1)
            % 8.5. it seems that the randomness wasn't enough
             R_i{1}(ii,:)=R_i{1}(ii,:)+0.5*Input.p(1)*rand([1 1]);
%              R_i{1}(ii,:)=R_i{1}(ii,:)+*Input.p(1)*rand([1 1]);
        end
        end
%         R_i{1}=R_i{1}+Input.p(1)*rand(size(R_i{1}));
            
%         save tmp_1D
        % build Nx2 ri as output position vector
        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        
%         figure; scatter(Array.r_i(:,1),Array.r_i(:,2)); axis square
        
        switch Input.sample.type
            case 'stiched_array'
                % shift the quadrants of the array to mimick stitching
                % error:                
                Ix=find(Array.r_i(:,1)>array_max_xlength/2);    % x sticthing:
                Iy=find(Array.r_i(:,2)>array_max_ylength/2);    % y sticthing:
                Array.r_i(Ix,1)=Input.sample.sticthing(1)+Array.r_i(Ix,1);
                Array.r_i(Ix,2)=Input.sample.sticthing(2)+Array.r_i(Ix,2);             
                Array.r_i(Iy,2)=Input.sample.sticthing(2)+Array.r_i(Iy,2);
        
            otherwise
        end
               
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            %max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end 
    case 'PC_mirrors'      
        %not supported for NPs in z direction
        if length(Input.arrayN) >2
            mirror_thickness = Input.PC_mirror_thickness;
            p = Input.p(1);
            N = Input.arrayN(1);
            wall_p = Input.PC_mirror_p(1);
            offset = Input.PC_mirror_offset(1); 

            n_mid = N - 2*mirror_thickness - 1;
            z1 = (0:mirror_thickness) .* wall_p;           
            z2 = z1(end) + offset + [1:n_mid] .* p;            
            z3 = z2(end) + offset + (1:mirror_thickness) .* wall_p;
            r_xi = [z1 z2 z3];

            p = Input.p(2);
            N = Input.arrayN(2);
            wall_p = Input.PC_mirror_p(2);
            offset = Input.PC_mirror_offset(2); 

            n_mid = N - 2*mirror_thickness - 1;
            z1 = (0:mirror_thickness) .* wall_p;           
            z2 = z1(end) + offset + [1:n_mid] .* p;            
            z3 = z2(end) + offset + (1:mirror_thickness) .* wall_p;
            r_yi = [z1 z2 z3];

            r_zi = 0;
       
            for ll = 1:(Input.array_shifts_y)
                
                 r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)=r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)+Input.p(2)/2;
            end
            for ll = 1:(Input.array_shifts_x)
                
                 r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)=r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)+Input.p(2)/2;
            end

            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
            
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;       
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end              
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);

        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        

               
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            %max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end 
    case 'PC_mirror_reflection_test'      
        %not supported for NPs in z direction
        % Creates a PC mirror in the middle of the array in y. Use only for
        % if N(2)-mirror_thickess is even.
        mirror_thickness = Input.PC_mirror_thickness;
        p = Input.p(1);
        N = Input.arrayN(1);
        wall_p = Input.PC_mirror_p(1);
        offset = Input.PC_mirror_offset(1); 

        n_1 = (N - mirror_thickness) / 2;
     
        z1 = [0:ceil(n_1-1)] * p;
        z2 = z1(end) + offset + p + [0:(mirror_thickness-1)] * wall_p;
        z3 = z2(end) + offset + [1:floor(n_1)] * p;
        r_xi = [z1 z2 z3];

        p = Input.p(2);
        N = Input.arrayN(2);
        wall_p = Input.PC_mirror_p(2);
        offset = Input.PC_mirror_offset(2); 

        n_1 = (N - mirror_thickness) / 2;
     
      

        z1 = [0:ceil(n_1-1)] * p;
        z2 = z1(end) + offset + p + [0:(mirror_thickness-1)] * wall_p;
        z3 = z2(end) + offset + [1:floor(n_1)] * p;
        r_yi = [z1 z2 z3];

        r_zi = 0;
   
        for ll = 1:(Input.array_shifts_y)
            
             r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)=r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)+Input.p(2)/2;
        end
        for ll = 1:(Input.array_shifts_x)
            
             r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)=r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)+Input.p(2)/2;
        end

        Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
        
        
        
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);

        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        

               
        %% generalize dipole periods to 3D case, if not already done so:
        %
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}] = meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end 
        %}
    case 'PC_mirror_reflection_test_reference'     
        % USE ONLY WITH PC MIRRROR REFLECTION TEST
        %not supported for NPs in z direction
        % Creates a PC mirror in the middle of the array in y. Use only for
        % if N(2)-mirror_thickess is even.
        mirror_thickness = Input.PC_mirror_thickness;
        p = Input.p(1);
        N = Input.arrayN(1);
        wall_p = Input.PC_mirror_p(1);
        offset = Input.PC_mirror_offset(1); 

        n_1 = (N - mirror_thickness) / 2;
     
        z1 = [0:ceil(n_1-1)] * p;
        z2 = z1(end) + offset + p + [0:(mirror_thickness-1)] * wall_p;
        z3 = z2(end) + offset + [1:floor(n_1)] * p;

        r_xi = [z1 z2 z3];

        p = Input.p(2);
        N = Input.arrayN(2);
        wall_p = Input.PC_mirror_p(2);
        offset = Input.PC_mirror_offset(2); 

        n_1 = (N - mirror_thickness) / 2;
     
        z1 = [0:ceil(n_1-1)] * p;
        z2 = z1(end) + offset + p + [0:(mirror_thickness-1)] * wall_p;
        z3 = z2(end) + offset + [1:floor(n_1)] * p;
        r_yi = [z1 z3];

        r_zi = 0;
   
        for ll = 1:(Input.array_shifts_y)
            
             r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)=r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)+Input.p(2)/2;
        end
        for ll = 1:(Input.array_shifts_x)
            
             r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)=r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)+Input.p(2)/2;
        end

        Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
        
         
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);

        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        

               
        %% generalize dipole periods to 3D case, if not already done so:
%
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end 
        %}
    case 'PC_middle'  
        % NOT IMPLEMENTED
        error('PC_middle array type is not implemented')
        %not supported for NPs in z direction
        if length(Input.arrayN) >2
            mirror_thickness = Input.PC_mirror_thickness;
            p = Input.p(1);
            N = Input.arrayN(1);
            wall_p = Input.PC_mirror_p(2);
            offset = Input.PC_mirror_offset(1); 

            n_mid = N - 2*mirror_thickness - 1;
            z1 = (0:mirror_thickness) .* wall_p;           
            z2 = z1(end) + offset + [1:n_mid] .* p;            
            z3 = z2(end) + offset + (1:mirror_thickness) .* wall_p;
            r_xi = [z1 z2 z3];

            p = Input.p(2);
            N = Input.arrayN(2);
            wall_p = Input.PC_mirror_p(2);
            offset = Input.PC_mirror_offset(2); 

            n_mid = N - 2*mirror_thickness - 1;
            z1 = (0:mirror_thickness) .* wall_p;           
            z2 = z1(end) + offset + [1:n_mid] .* p;            
            z3 = z2(end) + offset + (1:mirror_thickness) .* wall_p;
            r_yi = [z1 z2 z3];

            r_zi = 0;
       
            for ll = 1:(Input.array_shifts_y)
                
                 r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)=r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)+Input.p(2)/2;
            end
            for ll = 1:(Input.array_shifts_x)
                
                 r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)=r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)+Input.p(2)/2;
            end

            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
            
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;       
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end              
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);

        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        

               
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            %max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end 
    case 'double_rectangular'
        if length(Input.arrayN) >2
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
            array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
            array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);            
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
            r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
            r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
            
            rxi2 = Input.second_array_offset(1); 
            ryi2 = Input.second_array_offset(2); 
            rzi2 = Input.second_array_offset(3); 

            for ll = 1:(Input.array_shifts_y)
                
                 r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)=r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)+Input.p(2)/2;
            end
            for ll = 1:(Input.array_shifts_x)
                
                 r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)=r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)+Input.p(2)/2;
            end

            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
            
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;       
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end              
        [R_i1{1},R_i1{2},R_i1{3}] = meshgrid(r_xi,r_yi,r_zi);
        [R_i2{1},R_i2{2},R_i2{3}] = meshgrid(r_xi+rxi2,r_yi+ryi2,r_zi+rzi2);  

        %% add slight randomness to the the other (x) direction coordinates here:     
        if 0
        for ii=1:Input.Nmax(1)
            % 8.5. it seems that the randomness wasn't enough
             R_i{1}(ii,:)=R_i{1}(ii,:)+0.5*Input.p(1)*rand([1 1]);
%              R_i{1}(ii,:)=R_i{1}(ii,:)+*Input.p(1)*rand([1 1]);
        end
        end
%         R_i{1}=R_i{1}+Input.p(1)*rand(size(R_i{1}));
            
%         save tmp_1D
        % build Nx2 ri as output position vector
        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
            R_i{ii} = [R_i1{ii}, R_i2{ii}];
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        
%         figure; scatter(Array.r_i(:,1),Array.r_i(:,2)); axis square
        
        switch Input.sample.type
            case 'stiched_array'
                % shift the quadrants of the array to mimick stitching
                % error:                
                Ix=find(Array.r_i(:,1)>array_max_xlength/2);    % x sticthing:
                Iy=find(Array.r_i(:,2)>array_max_ylength/2);    % y sticthing:
                Array.r_i(Ix,1)=Input.sample.sticthing(1)+Array.r_i(Ix,1);
                Array.r_i(Ix,2)=Input.sample.sticthing(2)+Array.r_i(Ix,2);             
                Array.r_i(Iy,2)=Input.sample.sticthing(2)+Array.r_i(Iy,2);
        
            otherwise
        end
               
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            %max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end 

    case {'Outerwall', "HG_wall", "Inverted_Gaussian","Gaussian"}      
        
        % maximum array dimensions
        array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
        array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
        array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);            
        
        % x-, y- and z-coordinates
        r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
        r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
        r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
        
        for ll = 1:(Input.array_shifts_y)
            
             r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)=r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)+Input.p(2)/2;
        end
        for ll = 1:(Input.array_shifts_x)
            
             r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)=r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)+Input.p(2)/2;
        end

        Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
            
        
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
              

        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        
%         figure; scatter(Array.r_i(:,1),Array.r_i(:,2)); axis square

               
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            %max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end  
    case '1Dsquare_halfshifts2' 
        if length(Input.arrayN) >2
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
            array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
            array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);            
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
            r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
            r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
            
            for ll = 1:(Input.array_shifts_y)
                
                 r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)=r_yi(ceil(length(r_yi)/(Input.array_shifts_y+1)*ll):end)-Input.p(2)/2;
            end
            for ll = 1:(Input.array_shifts_x)
                
                 r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)=r_xi(ceil(length(r_xi)/(Input.array_shifts_x+1)*ll):end)-Input.p(2)/2;
            end

            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
            
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;       
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end              
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
              
        %% add slight randomness to the the other (x) direction coordinates here:     
        if 0
        for ii=1:Input.Nmax(1)
            % 8.5. it seems that the randomness wasn't enough
             R_i{1}(ii,:)=R_i{1}(ii,:)+0.5*Input.p(1)*rand([1 1]);
%              R_i{1}(ii,:)=R_i{1}(ii,:)+*Input.p(1)*rand([1 1]);
        end
        end
%         R_i{1}=R_i{1}+Input.p(1)*rand(size(R_i{1}));
            
%         save tmp_1D
        % build Nx2 ri as output position vector
        for ii=1:3 % 3
            % Array.r_i(:,ii)=single(R_i{ii}(1:end)); 
             posed = pagetranspose(R_i{ii});
                Array.r_i(:,ii)=single(posed(1:end)); %transpose fixes dimension problems later, now it stays [x,y,z] when reshaped later 
        end
        
%         figure; scatter(Array.r_i(:,1),Array.r_i(:,2)); axis square
        
        switch Input.sample.type
            case 'stiched_array'
                % shift the quadrants of the array to mimick stitching
                % error:                
                Ix=find(Array.r_i(:,1)>array_max_xlength/2);    % x sticthing:
                Iy=find(Array.r_i(:,2)>array_max_ylength/2);    % y sticthing:
                Array.r_i(Ix,1)=Input.sample.sticthing(1)+Array.r_i(Ix,1);
                Array.r_i(Ix,2)=Input.sample.sticthing(2)+Array.r_i(Ix,2);             
                Array.r_i(Iy,2)=Input.sample.sticthing(2)+Array.r_i(Iy,2);
        
            otherwise
        end
               
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
                
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            %max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            max_size(ii)=Input.p(1,ii)*(Input.Nmax(ii)-1);
            
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            %Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            posed = pagetranspose(R_i{ii});
            Array.r_i_conv(:,ii)=single(posed(1:end));
        end       

        %         %% plot the array for testing:
        %         figure
        %         scatter(Array.r_i_conv(:,1),Array.r_i_conv(:,2),'x')
        %         set(gcf,'position',[139   156  1594   923])
        


    case 'fcc_square'
        
        %% Make checkerboad pattern:
        num_zero = 1e-3;    % numerical zero
%         num_zero = 1e-5;    % numerical zero, 21.8. trial
        nx=ceil(Input.arrayN/2);
        C_tmp=repmat([1 num_zero ; num_zero 1],[nx nx]);

        %% remove the additional line if ny is odd, since the eye(2) always makes a even array:
        if mod(Input.arrayN,2) == 0
            Array.checkerboard=C_tmp;
        else
            Array.checkerboard=C_tmp(1:(end-1),1:(end-1));
        end
        
        if length(Input.arrayN) >2
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
            array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
            array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
            r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
            r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
            
            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
        
        % build Nx2 ri as output position vector
        for ii=1:3 % 3
            Array.r_i(:,ii)=single(R_i{ii}(1:end));
        end
        
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
        
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
        end

    case 'square_random'
            
            Input.array_size=200e-6;           % physical array size
   
            [r_xi,r_yi]=Random_NP_sample_design_DDA(Input);
            [R_i{1},R_i{2}]=Random_NP_sample_design_DDA(Input);
            R_i{3}=zeros(size(R_i{1}));
            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(1) 1];
            
            % build Nx2 ri as output position vector
            for ii=1:3 % 3
                Array.r_i(:,ii)=single(R_i{ii}(1:end));
            end
            % build Nx3 ri as output position vector
            for ii=1:length(R_i)
                Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            end
 
    case 'square_perturbed'
        %% Payman's stuff here
        if length(Input.arrayN) >2
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
            array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
            array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN(1));
            r_yi=linspace(0,array_max_ylength,Input.arrayN(2));
            r_zi=linspace(0,array_max_zlength,Input.arrayN(3));
            
            Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
        else
            %% old codes:
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            % x-, y- and z-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            r_zi=0;
            
            Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
            
        end
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
        
        % build Nx2 ri as output position vector
        for ii=1:3 % 3
            %             Array.r_i(:,ii)=double(R_i{ii}(1:end));
            Array.r_i(:,ii)=single(R_i{ii}(1:end));
        end
        
        %% generalize dipole periods to 3D case, if not already done so:
        if length(Input.p) < 3
            Input.p(1,1,3)=Input.p(1,1,2);
        end
        
        for ii=1:length(Input.p) % depending on the dimensionality of the array
            max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
        end
        
    case 'hexagonal'
        % Make the hexagonal array from one rectangular array with px ~ py, 
        % and having scatting particles in a checkerboard pattern in the array:
        % here array line separation in y-direction is (equilateral triangular
        % lattice) h=sqrt(3)/2*r, r=r_xi(2)-r_xi(1), x-period is still px:
                
        % periodicities for the rectangular array:
        p_sq(1)=Input.p(1)/2;              % p_sq = 2*h = sqrt(3)*r, and px = py be definition
        p_sq(2)=sqrt(3)*Input.p(1)/2;      % p_sq = 2*h = sqrt(3)*r
        
        %% redefine the lattice size for hexagonal arrays, since we're superimposing two rectangular lattices:
        fun=@find_root_for_correct_hex_array;
%         Input.arrayN=floor(fzero(fun,Input.arrayN/2));
        Input.arrayN=ceil(fzero(fun,Input.arrayN/2));        
        [Array.checkerboard,Array.checkb_Innz]=Make_checkerboard_array(Input.arrayN,Input.arrayN);
              
        Input.arrayNx=size(Array.checkerboard,2);
        Input.arrayNy=size(Array.checkerboard,1);
        Input.arrayN=Input.arrayNx;
        
        % maximum array dimensions
        array_max_xlength=p_sq(1)*(Input.arrayNx-1);
        array_max_ylength=p_sq(2)*(Input.arrayNy-1);
             
        
        %% x-coordinates
        r_xi=linspace(0,array_max_xlength,Input.arrayNx);
        
        %% y-coordinates
        r_yi=linspace(0,array_max_ylength,Input.arrayNy);
        
        %% z-coordinates
        r_zi=0;
        
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
           
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            Array.r_i(:,ii)=double(R_i{ii}(1:end));
        end
        
        %% generalize dipole periods to 3D case, if not already done so:
        % but here with hexagonal arrays, px ~= py ! 
        p_sq(1)=Input.p(1)/2;              % p_sq = 2*h = sqrt(3)*r, and px = py be definition
        p_sq(2)=sqrt(3)*Input.p(1)/2;      % p_sq = 2*h = sqrt(3)*r
        p_sq(3)=p_sq(1);
            
        %% meshgridding should be easier and faster:
        Input.Nmax(1:3)=[Input.arrayNx Input.arrayNy 1];
        
        for ii=1:length(p_sq) % depending on the dimensionality of the array
            max_size(ii)=p_sq(ii)*(Input.Nmax(ii)-1);
            r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
        end
        
        clear R* r_*
    case 'stl_file'      
        
        %% load an .stl file
        % use a standard folder to keep the stl files:
        filename=[Input.stl_folder Input.stl_filename];
        [Array]=Load_stl_file(filename,Input.dipole_spacing);       % loads an .stl file
        
        %% meshgridding should be easier and faster:
        Input.Nmax(1:3)=Array.arrayN_ii;
               
        % build Nx3 ri_conv to be used in calculating Ax
        % basically it's r_i extended in [-MAX,MAX] size instead of [0,MAX]:
        r_i_cell=cell(3,1);
        for ii=1:3               
            r_i_cell{ii}=linspace(-Array.MAX_of_r_i(ii),Array.MAX_of_r_i(ii),2*Array.arrayN_ii(ii)-1);
        end
        
        %% form Arrays
        [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
        
        % build Nx3 ri as output position vector
        for ii=1:length(R_i)
            Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
        end
          
end
;
end

function [Oligomer,Input]=Make_3D_Oligomer_sub(Input)

switch Input.sample.oligomertype
    case 'yagi'
%          save yagi_array_tmp
         
         %% old codes:
         % maximum array dimensions
         array_max_xlength=Input.p(1)*(Input.arrayN-1);
         array_max_ylength=Input.p(2)*(Input.arrayN-1);
         
         % x-, y- and z-coordinates
         r_xi=linspace(0,array_max_xlength,Input.arrayN);
         r_yi=linspace(0,array_max_ylength,Input.arrayN);
         r_zi=0;
         
         Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];
         
         
         % assuming equidistant particle separation, r = px/5
         r=Input.Px(1)/5;
         
         %% jth particles, -j*r from original
         for jj=1:4
             r_xi_nu=r_xi - jj*r ;
             [R_i{1*jj},R_i{4+jj}]=meshgrid(r_xi_nu,r_yi);                     
         end
         R_nu{1}=[R_i{1} R_i{2} R_i{3} R_i{4}];
         R_nu{2}=[R_i{5} R_i{6} R_i{7} R_i{8}];
          
         % build Nx2 ri as output position vector
         for ii=1:3 % 3
             if ii==3
                 Oligomer.r_i(:,ii)=0;
             else
                 Oligomer.r_i(:,ii)=single(R_nu{ii}(1:end));
             end
         end
         
    case 'trimer'
        %% parameters
        N=3;     % number of scatterers making the oligomer
        Max_dim=5e-6;       % maximum size of the oligomer, should be smaller than possible periodicity of the array
        
        %% x-coordinates
        % assumed now to be spherically symmetric, i.e. coordinates on a sphere
        %   (x-x0)^2 + (y-y0)^2 = r^2 or in parameterized form
        %
        %  x = x0 + r*cos(theta)*sin(varphi)
        %  y = y0 + r*sin(theta)*sin(varphi) , where 0 < theta < 2*pi  and  0 < varphi < pi)
        %
        % and x0 and y0 describe the center of the lattice point
        % we could generalize to 3D with z and varphi to model e.g. helices
        
        x0=0;       % these should be compatible for replacing with Array structs array points
        y0=0;
        
        theta=deg2rad([0 120 240]);
        varphi=pi/2;
        r=Input.Px(1);
        
        r_xi=x0 + r*cos(theta).*sin(varphi);
        r_yi=y0 + r*sin(theta).*sin(varphi);
        
        %% Scatterer orientations
        theta_offset=deg2rad(21);
        Oligomer.theta=theta_offset + deg2rad([0 120 240]);
        
        phi_offset=deg2rad(15);
        Oligomer.phi=phi_offset + deg2rad([0 0 0]);
        
        Oligomer.psi=Oligomer.phi;
        Oligomer.psi(:)=0;
        
        %     function r_xi_fun(x0,theta,varphi)
        %         % change the sizes of x0 and theta to be compatible with each other
        %
        %         repmat(x0,[1 size(theta)]) + repmat(r*cos(theta).*sin(varphi),[size(x0,1) size(x0,2)]));
        %
        %
        %     end
        % r_xi_fun(x0,theta,varphi)
        
        
        %%
        % [R_i{1},R_i{2}]=meshgrid(r_xi,r_yi);
        R_i{1} = r_xi.';
        R_i{2} = r_yi.';
        
        % build Nx2 ri as output position vector
        for ii=1:2 % 3
            Oligomer.r_i(:,ii)=single(R_i{ii}(1:end));
        end
        %%
    case 'line'
        %% parameters
        if Input.arrayN > 10
                   
            % maximum array dimensions
            array_max_xlength=Input.p(1)*(Input.arrayN-1);
            array_max_ylength=Input.p(2)*(Input.arrayN-1);
            
            %% x-coordinates
            r_xi=linspace(0,array_max_xlength,Input.arrayN);
            
            % depending on the oddness/evenness of N, offset the array
            % %0...max  - max/2 = -max/2..max/2 BUT need to force particle to
            % zero point
            if mod(Input.arrayN,2) == 0
                r_xi=r_xi-array_max_xlength/2-Input.p(1)/2;
            else
                r_xi=r_xi-array_max_xlength/2;
            end
            
            %% y-coordinates
            r_yi=linspace(0,array_max_ylength,Input.arrayN);
            
            if mod(Input.arrayN,2) == 0
                r_yi=r_yi-array_max_ylength/2-Input.p(2)/2;
            else
                r_yi=r_yi-array_max_ylength/2;
            end
            
            %%
            [R_i{1},R_i{2}]=meshgrid(r_xi,r_yi);
            
            %% form now N times new arrays where the locations are slightly shifted
            N=3;                % number of scatterers making the oligomer
            Max_dim=5e-6;       % maximum size of the oligomer, should be smaller than possible periodicity of the array
            
            % form single oligomer relative coordinates:
%             r=Input.Px(1)*0.1;  
            r=Input.Px(1);
            
            %% 2. particles, -r from originals
            r_xi_nu=r_xi -r ;                       
            [R_i{3},R_i{4}]=meshgrid(r_xi_nu,r_yi);
            
            %%
            %% 3. particles, +r from originals
            r_xi_nu=r_xi +r ;           
            [R_i{5},R_i{6}]=meshgrid(r_xi_nu,r_yi);
              
            R_nu{1}=[R_i{1} R_i{3} R_i{5}];
            
            R_nu{2}=[R_i{2} R_i{4} R_i{6}];
            
%             %% z-coordinates
%             r_zi=0;
        
            %%
%             [R_i{1},R_i{2},R_i{3}]=meshgrid(r_xi,r_yi,r_zi);
            % build Nx2 ri as output position vector
            for ii=1:3 % 3
                if ii==3
                    Oligomer.r_i(:,ii)=0;
                else
                    Oligomer.r_i(:,ii)=single(R_nu{ii}(1:end));
                end
            end
            
            %%
        elseif Input.arrayN == 3
            N=Input.arrayN;                % number of scatterers making the oligomer
            Max_dim=5e-6;       % maximum size of the oligomer, should be smaller than possible periodicity of the array
                        
            x0=0;       % these should be compatible for replacing with Array structs array points
            y0=0;
            
            r=Input.Px(1);
            
            r_xi=x0 + [-r 0 +r];
%             r_xi=x0 + [-2*r -r 0 +r +2*r];
            r_yi=repmat(y0,size(r_xi));
            
            
            %%
            % [R_i{1},R_i{2}]=meshgrid(r_xi,r_yi);
            R_i{1} = r_xi;
            R_i{2} = r_yi;
            
            % build Nx3 ri as output position vector
            for ii=1:3 % 3
                if ii==3
                    Oligomer.r_i(:,ii)=0;
                else
                    Oligomer.r_i(:,ii)=single(R_i{ii}(1:end));
                end
            end
            
            Oligomer.r_i

            %% generalize dipole periods to 3D case, if not already done so:
            if length(Input.p) < 3
                Input.p(1,1,3)=Input.p(1,1,2);
            end
            
            %% meshgridding should be easier and faster:
            Input.Nmax(1:3)=[Input.arrayN 1 1];
            
            for ii=1:length(Input.p) % depending on the dimensionality of the array
                max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
                r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
            end
            
            %% form Arrays
            [R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});
            
            % build Nx3 ri as output position vector
            for ii=1:length(R_i)
                Oligomer.r_i_conv(:,ii)=single(R_i{ii}(1:end));
            end
            
              elseif Input.arrayN == 5
            N=5;                % number of scatterers making the oligomer
            Max_dim=5e-6;       % maximum size of the oligomer, should be smaller than possible periodicity of the array
            
            %% x-coordinates
            % assumed now to be spherically symmetric, i.e. coordinates on a sphere
            %   (x-x0)^2 + (y-y0)^2 = r^2 or in parameterized form
            %
            %  x = x0 + r*cos(theta)*sin(varphi)
            %  y = y0 + r*sin(theta)*sin(varphi) , where 0 < theta < 2*pi  and  0 < varphi < pi)
            %
            % and x0 and y0 describe the center of the lattice point
            % we could generalize to 3D with z and varphi to model e.g. helices
            
            x0=0;       % these should be compatible for replacing with Array structs array points
            y0=0;
            
            r=Input.Px(1)*0.1;
            
%             r_xi=x0 + [-r 0 +r];
            r_xi=x0 + [-2*r -r 0 +r +2*r];
            r_yi=repmat(y0,size(r_xi));
            
            
            %%
            % [R_i{1},R_i{2}]=meshgrid(r_xi,r_yi);
            R_i{1} = r_xi;
            R_i{2} = r_yi;
            
            % build Nx2 ri as output position vector
%             for ii=1:2 % 3
%                 Oligomer.r_i(:,ii)=single(R_i{ii}(1:end));
%             end
            for ii=1:3 % 3
                if ii==3
                    Oligomer.r_i(:,ii)=0;
                else
                    Oligomer.r_i(:,ii)=single(R_i{ii}(1:end));
                end
            end
            %%
        end
end

end

function [Array,Input]=Make_Random_Array_sub(Input);   

% save tmp_ran_array

if length(Input.arrayN) >2
    % maximum array dimensions
    array_max_xlength=Input.p(1)*(Input.arrayN(1)-1);
    array_max_ylength=Input.p(2)*(Input.arrayN(2)-1);
    array_max_zlength=Input.p(3)*(Input.arrayN(3)-1);
    
    disp('3D array')
    %     % x-, y- and z-coordinates
    %     r_xi=linspace(0,array_max_xlength,Input.arrayN);
    %     r_yi=linspace(0,array_max_ylength,Input.arrayN);
    %     r_zi=0;
    
    %% Random surface:
    rng(Input.ran_mod*50);
    R_i{1}=array_max_xlength*rand(Input.arrayN(1),Input.arrayN(2),Input.arrayN(3));
    R_i{2}=array_max_ylength*rand(Input.arrayN(1),Input.arrayN(2),Input.arrayN(3));
    R_i{3}=array_max_zlength*rand(Input.arrayN(1),Input.arrayN(2),Input.arrayN(3));
        
    Input.Nmax(1:3)=[Input.arrayN(1) Input.arrayN(2) Input.arrayN(3)];
else
    disp('2D array')
    %% old codes:
    % maximum array dimensions
    array_max_xlength=Input.p(1)*(Input.arrayN-1);
    array_max_ylength=Input.p(2)*(Input.arrayN-1);
    
%     % x-, y- and z-coordinates
%     r_xi=linspace(0,array_max_xlength,Input.arrayN);
%     r_yi=linspace(0,array_max_ylength,Input.arrayN);
%     r_zi=0;
    
    %% Random surface:
    rng(50);
    
%     pm1 = rand([Input.arrayN Input.arrayN]);
%     pm1(pm1<0.5) = -1;
%     pm1(pm1>=0.5) = 1;
% 
%     pm2 = rand([Input.arrayN Input.arrayN]);
%     pm2(pm2<0.5) = -1;
%     pm2(pm2>=0.5) = 1;
%     
%     R_i{1}=pm1.*0.5*array_max_xlength*rand(Input.arrayN,Input.arrayN);
%     R_i{2}=pm2.*0.5*array_max_ylength*rand(Input.arrayN,Input.arrayN);
    
    R_i{1}=array_max_xlength*rand(Input.arrayN,Input.arrayN);
    R_i{2}=array_max_ylength*rand(Input.arrayN,Input.arrayN);
    
    R_i{3}=zeros(Input.arrayN,Input.arrayN);
    
    Input.Nmax(1:3)=[Input.arrayN Input.arrayN 1];    
end

% build Nx2 ri as output position vector
for ii=1:3 % 3
    Array.r_i(:,ii)=single(R_i{ii}(1:end));
end

% save scatter_tmp
% sketch the beam into the scatter plot:
% r=3e-6;
% ang=0:0.01:2*pi; 
% xp=r*cos(ang);
% yp=r*sin(ang);

% save sample
figure
if length(Input.arrayN) >2    
    scatter3(Array.r_i(:,1),Array.r_i(:,2),Array.r_i(:,3))
else
    
    scatter(Array.r_i(:,1),Array.r_i(:,2))
    hold on
    plot((max(Array.r_i(:,1))/2)+xp,(max(Array.r_i(:,2))/2)+yp);
end

%% generalize dipole periods to 3D case, if not already done so:
if length(Input.p) < 3
    Input.p(1,1,3)=Input.p(1,1,2);
end

for ii=1:length(Input.p) % depending on the dimensionality of the array
    max_size(ii)=Input.p(1,1,ii)*(Input.Nmax(ii)-1);
    r_i_cell{ii}=linspace(-max_size(ii),max_size(ii),2*Input.Nmax(ii)-1);
end

%% form Arrays
[R_i{1},R_i{2},R_i{3}]=meshgrid(r_i_cell{1},r_i_cell{2},r_i_cell{3});

% build Nx3 ri as output position vector
for ii=1:length(R_i)
    Array.r_i_conv(:,ii)=single(R_i{ii}(1:end));
end

end