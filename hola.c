    if(strcmp(window_function, "yes") == 0 && strcmp(shuffle,"no") == 0){
    
        if(strcmp(type_of_code, "rustico") == 0){sprintf(name_wink_out,"%s/Window_%s.txt",name_path_out,name_id);}
    if(strcmp(type_of_code, "rusticoX") == 0){
    sprintf(name_wink_out,"%s/WindowAA_%s.txt",name_path_out,name_id);
    sprintf(name_winkBB_out,"%s/WindowBB_%s.txt",name_path_out,name_id);
    sprintf(name_winkAB_out,"%s/WindowAB_%s.txt",name_path_out,name_id);
    }
        

    f=fopen(name_wink_out,"w");
    if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_wink_out);return 0;}
    fclose(f);
if(strcmp(type_of_code, "rusticoX") == 0){
    f=fopen(name_winkBB_out,"w");
    if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_wink_out);return 0;}
    fclose(f);

    f=fopen(name_winkAB_out,"w");
    if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_wink_out);return 0;}
    fclose(f);
}

    printf("Computing the RR counts using %d processors. This can take a while...\n",n_lines_parallel);
//window_mask_function_RRcount(name_wink_out,name_randoms_in,window_norm_bin,deltaS_window,percentage_randoms_window,yamamoto4window,parameter_value,header,L2-L1,n_lines_parallel );

 window_mask_function_RRcount(name_wink_out,name_winkBB_out,name_winkAB_out,name_randoms_in,name_randomsB_in,window_norm_bin,deltaS_window,percentage_randoms_window,yamamoto4window,parameter_value,parameter_valueB,header,L2-L1,n_lines_parallel,type_of_code);


        if(strcmp(type_of_code, "rustico") == 0){
    if(Ndata==0){printf("Window function counts completed. No data file provided for power spectrum computation. Exiting now...\n");return 0;}
    if(Ndata>0){printf("Window function counts completed.\n\n");}
            parameter_value[3]=Ndata;
        }
        if(strcmp(type_of_code, "rusticoX") == 0){
            if(Ndata==0 || NdataB==0){printf("Window function counts completed. No data file provided for power spectrum computation. Exiting now...\n");return 0;}
            if(Ndata>0 && NdataB>0){printf("Window function counts completed.\n\n");}
            parameter_value[3]=Ndata;
            parameter_valueB[3]=NdataB;

        }

    }

