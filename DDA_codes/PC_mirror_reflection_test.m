function PC_mirror_reflection_test(Input,Array)
        %% first calculations with the mirror in place
        
        [P1,C1,Field,Array_out] = Calculate_3D_DDA_generic(Input,Array);
        %% remove the mirror for reference
        
        Input.array_type = 'PC_mirror_reflection_test_reference';
        [Array,Input]=Make_generic_Array(Input);
        Input.arrayN = [Input.arrayN(1) Input.arrayN(2)-Input.PC_mirror_thickness-1 1];

        [Array]=NLO_DDA_Make_3D_polarizability(Array,Input);

          [P2,C2,Field,Array_out2] = Calculate_3D_DDA_generic(Input,Array);

         %% Calculate the difference between reference and with mirror
Array_out2

end