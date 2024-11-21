% This function is written by Mingliang Wu (email: wuml_neu@163.com)

function [ps,pf,mertic_iter] = MOCDA(func_name,LB,UB,dim,N,Max_fevs,n_var,repoint)
rand('state', sum(100*clock)); 
load (strcat([func_name,'_Reference_PSPF_data']));
%% Initialization parameters
active_level = 1.0; aa = 0.25; bb = 0.5; cc = 0.75; % Active_level is best left unchanged because after experimentation we discarded SCAE
FES = 0; 
DD = 1; AAA = 5; BBB = 1;
ind.Position=[]; 
ind.Cost=[];
ind.Rank=[];
ind.DominationSet=[];
ind.DominatedCount=[];
ind.CrowdingDistance=[];
pop=repmat(ind,N,1);
for ii=1:N % Initializing populations
    pop(ii).Position=LB+(UB-LB).*rand(1,dim);
    Population_decs(ii,:) = pop(ii).Position;
end
for ii=1:N % Calculate the fitness
    aaabbb = feval(func_name,pop(ii).Position);
    Population_objs(ii,:) = aaabbb;
end  
    Parentstar = [Population_decs,Population_objs];
    K = 10;
    Parentstar = non_domination_scd_kmeans_sort(Parentstar,n_var,dim,K); %Cluster
for i1 = 1:size(Parentstar,1) % Reassignment
    pop(i1).Position = Parentstar(i1,1:dim);
    pop(i1).Cost = Parentstar(i1,dim+1:dim+n_var)'; 
    pop(i1).Rank = Parentstar(i1,dim+n_var+1); 
    pop(i1).CrowdingDistance = Parentstar(i1,dim+n_var+3); 
end
%% Loop
i = 1; % Iteration counter
while FES<Max_fevs 
    pop_old = pop;
    newpop=repmat(ind,N,1);
    for ii=1:N
        P(ii,:) = pop(ii).Position;
    end
    F = 1:1:N;
    for iy = 1:size(pop,1)
        lp(iy) = pop(iy).Rank;
    end
    for i5 = 1:max(lp)
        eval(['s = ',num2str(i5),';'])
        eval(['[~,l] = find(lp==s); CD',num2str(i5),' = P(l,:); AF = F(l); [~,FM',num2str(i5),'] = sort(AF);'])
    end
    NP = zeros(N,dim); 
    for ii = 1:N
        inter = P(ii,:); 
        if rand > active_level
            if lp(ii) == 1; NP(ii,:) = inter; end % This step doesn't change the individual, so it doesn't update the FES
            if lp(ii) > 1
                RR = rand;
                for i5 = 1:max(lp)
                    eval(['if  s == ',num2str(i5),'; CD1 = CD',num2str(i5),'; FM1 = FM',num2str(i5),'; end'])
                end   

                if RR<aa % M
                    MP = FM1(1); NP(ii,:) = DD*rand(1,dim).*(CD1(MP,:)-inter)+inter;
                    FES = FES + 1;
                end
                if  aa<=RR && RR<bb % E
                    guimo = numel(FM1);
                    interinter = CD1(FM1(randperm(ceil(AAA/50*guimo),ceil(BBB/50*guimo))),:);
                    if size(interinter,1) > 1
                        NP(ii,:) = inter + DD*rand(1,dim).*(sum(interinter) - ceil(BBB/50*guimo)*inter);
                    else
                        NP(ii,:) = inter + DD*rand(1,dim).*((interinter) - ceil(BBB/50*guimo)*inter);
                    end
                    FES = FES + 1;
                end   
                if bb<=RR && RR<cc % W
                    guimo = numel(FM1);
                    xulie = ceil(AAA/50*guimo)+1;
                    if max(xulie) > numel(FM1)
                        xulie = xulie - (max(xulie) - FM1);
                    end 
                    MP = FM1(xulie); NP(ii,:) = DD*rand(1,dim).*(CD1(MP,:)-inter)+inter;
                    FES = FES + 1;
                end
                if cc<=RR && RR<=1 % U
                    guimo = numel(FM1);
                    xulie = randperm(ceil(AAA/50*guimo),ceil(BBB/50*guimo))+ceil(AAA/50*guimo);
                    if max(xulie) > numel(FM1)
                        xulie = xulie - (max(xulie) - FM1);
                    end 
                    interinter = CD1(FM1(xulie),:);
                    if size(interinter,1) > 1
                        NP(ii,:) = inter + DD*rand(1,dim).*(sum(interinter) - ceil(BBB/50*guimo)*inter);
                    else
                        NP(ii,:) = inter + DD*rand(1,dim).*((interinter) - ceil(BBB/50*guimo)*inter);
                    end
                    FES = FES + 1;
                end
            end
        else 
             s = randi(max(lp)); % 3) in helpers selection mechanism
             for i5 = 1:max(lp) 
                 eval(['if  s == ',num2str(i5),'; CD = CD',num2str(i5),'; FM = FM',num2str(i5),'; end'])
             end
             eval(['guimo = numel(FM',num2str(s),');'])
             interinter = CD(FM(randperm(ceil(AAA/50*guimo),ceil(BBB/50*guimo))),:);
             if size(interinter,1) > 1
                 NP(ii,:) = inter + DD*rand(1,dim).*(sum(interinter) - ceil(BBB/50*guimo)*inter);
             else
                 NP(ii,:) = inter + DD*rand(1,dim).*((interinter) - ceil(BBB/50*guimo)*inter);
             end
             FES = FES + 1;
        end
    end
    for ioi = 1:N % Repair of infeasible solutions
        for joj = 1:dim
            NP(ioi,joj) = rand*(UB(1,joj)-LB(1,joj)) + LB(1,joj);
        end
    end
    for ii=1:N % Calculating the fitness of new population
        newpop(ii).Position = NP(ii,:);
        aaabbb = feval(func_name,NP(ii,:));
        if size(aaabbb,1) == 1
            aaabbb = aaabbb';
        end
        newpop(ii).Cost = aaabbb;
    end
    pop=[pop;newpop]; % Combining old and new populations
    for ii=1:size(pop,1)
        Population_decs(ii,:) = pop(ii).Position;
        Population_objs(ii,:) = pop(ii).Cost;
    end
    Parentstar = [Population_decs,Population_objs];
    K = 10;
    Parentstar = non_domination_scd_kmeans_sort(Parentstar,n_var,dim,K); %聚类
    for i1 = 1:size(Parentstar,1) % Reassignment
        pop(i1).Position = Parentstar(i1,1:dim);
        pop(i1).Cost = Parentstar(i1,dim+1:dim+n_var)'; 
        pop(i1).Rank = Parentstar(i1,dim+n_var+1); 
        pop(i1).CrowdingDistance = Parentstar(i1,dim+n_var+3); 
    end
    pop = pop(1:N); % Retain the first half of good individuals
%% Recording (in-loop)
    Rank_take = []; 
    for i3 = 1:size(pop,1)
        if pop(i3).Rank == 1
            Rank_take = [Rank_take,i3];
        end
    end
    pf_max_serial = numel(Rank_take); % Get the maximum of the pf

    for i4 = 1:pf_max_serial
        pf(i4,:) = pop(i4).Cost;
        ps(i4,:) = pop(i4).Position;
    end
    % If the PS, PF and repoint cannot be known in advance, you can omit the following steps, 
    % which will have little impact on the performance of the algorithm and even get better, 
    % but remember to change the script output.
    pf_cell{i} = pf;
    ps_cell{i} = ps;
    hyp = Hypervolume_calculation(pf,repoint); 
    IGDx = IGD_calculation(ps,PS);
    IGDf = IGD_calculation(pf,PF); 
    CR = CR_calculation(ps,PS); 
    PSP = CR/IGDx; 
    mertic_iter(i,:) = [1./PSP,1./hyp,IGDx,IGDf];
    if i > 100000 % Retain current generation evolution if overall progress has been made 
        if mertic_iter(i-1,2) > 0
            if sum(mertic_iter(i,:)) > sum(mertic_iter(i-1,:))
                mertic_iter(i,:) = mertic_iter(i-1,:);
                pop = pop_old;
                pf_cell{i} = pf_cell{i-1};
                ps_cell{i} = ps_cell{i-1};
            end
        end
    end
    i = i + 1; % Update Iteration Counter
end
%% Output
pf = pf_cell{end};
ps = ps_cell{end};

