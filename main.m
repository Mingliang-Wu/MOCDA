clc
clear
close
addpath(genpath('.\MM_testfunctions\'));
addpath(genpath('.\Indicator_calculation\'));
N_function = 24;% number of test function
runtimes = 21;  % odd number
    for i_func = 1:N_function 
        i_func;
        [fname,n_var,n_obj,xl,xu,repoint,N_ops] = switch_func(i_func);
        %% Load reference PS and PF data
        load  (strcat([fname,'_Reference_PSPF_data']));
        %% Initialize the population size and the maximum evaluations
        popsize = 200*N_ops;
        Max_fevs = 10000*N_ops;
        Max_Gen = fix(Max_fevs/popsize);
        Algorithm = Algorithm_Function();
        runtimes = 21; 
        for j=1:runtimes
            j
                eval(['[ps,pf,mertic_iter] = ',Algorithm{1},'(fname,xl,xu,n_obj,popsize,Max_fevs,n_var,repoint);']); 
                hyp = Hypervolume_calculation(pf,repoint); 
                IGDx = IGD_calculation(ps,PS); 
                IGDf = IGD_calculation(pf,PF); CR = CR_calculation(ps,PS); PSP = CR/IGDx; 
                Algorithm_name = 'MOCDA_improve'; 
                eval(['Indicator.',Algorithm_name,'(j,:) = [1./PSP,1./hyp,IGDx,IGDf]; PSdata.',Algorithm_name,'{j} = ps; PFdata.',Algorithm_name,'{j} = pf;']); 
                mertic_iter_cell{j} = mertic_iter; 
                clear ps pf hyp IGDx IGDf CR PSP 
        end
        AvgOfRankOfProblems_cell{i_func} = FriedmanTest (Indicator,Algorithm);
        save AvgOfRankOfProblems_cell AvgOfRankOfProblems_cell
        eval (['save results_',num2str(i_func),';']) 
    end