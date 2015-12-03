function [seq,population,plotterjuk] = EvolvePop(popfile,ngen, mutrate, recombrate, recombsize)

[header,pop] = fastaread(popfile);

bases={'A','C','G','T'};
basesint=[1 2 3 4];

popsize=length(header);
pop2=cell2mat(pop);
pop3=nt2int(pop2);
seqsize=length(pop3)/popsize;

seq=zeros(1, seqsize);

for z=1:seqsize
    seq(1,z)=pop3(z);
end

aux = 1;
j = 1;
for i = 1:length(pop3)

    if i > seqsize*j;
        j = j+1;
        aux = 1;
    end
    population(j,aux) = pop3(i);
    aux = aux + 1;
end

%% Evolution throughout generations
k=ngen;                                           % Number of iterations
identity = zeros(sum(1:popsize-1),1);             % Vector with values from Sequence Identity model on a given generation
jcantor  = zeros(sum(1:popsize-1),1);             % Vector with values from Jukes-Cantor model on a given generation
xplot    = ones(sum(1:popsize-1),1);              % Vector with values from the generations
triup    = triu(ones(popsize,popsize),1);
plotterseq  = zeros(ngen,1);
plotterjuk  = zeros(ngen,1);

for g=1:ngen,
    %% Mutations

    for m = 1:popsize
        % Evaluates if occurs a mutation
        rint = rand(1,1);
        if rint <= mutrate

            % Builds, randomly, the position in which a mutation happens

            posmut = randi([1,seqsize],1);

            % Builds a vector of the bases that can suffer a mutation
            newbases = basesint(basesint~=population(m,posmut));

            % Changes the nucleotide on the population matrix (posx,posy)
            population(m,posmut) = newbases(randi(size(newbases,2),1));

        end

    end

    %% Recombinations
    popperm = randperm(popsize);                   % Position permutation vector

    for r=1:popsize

        % Evaluates if occurs recombination
        rint = rand(1,1);
        if rint <= recombrate

            % Bulids, randomly, the position and individual in which occurs recombiantion

            indrec = randi([1,popsize],1);

            % Recombination only happens if final sequence is different from the initial
            if indrec ~= popperm(r)

                % Builds the position, on the initial sequence, where recombination occurs
                rposi=randi([1,seqsize-(recombsize-1)],1);

                % Builds the position, on the final sequence, where recombination occurs
                rposf=randi([1,seqsize-(recombsize-1)],1);

                % Vector with the sequece to recombine
                frag    = zeros(1,recombsize);

                for f    = 1:recombsize
                    frag = population(indrec,rposi:(rposi+recombsize-1));
                end

                % Recombines on the final sequence
                population(popperm(r),rposf:(rposf+recombsize-1)) = frag;
            end
        end
    end

    %% Builds distance using Sequence Identity (%) Model and  Jukes & Cantor Model
    identaux = zeros(popsize,popsize);
    jcantaux = zeros(popsize,popsize);

    for a = 1:popsize
        for b = a+1:popsize

            p = sum(population(a,:) ~= population(b,:))/seqsize;
            % Sequence Identity
            identaux(a,b) = p;
            % Jukes & Cantor
            jcantaux(a,b) = real(-(3/4)*log(1-(4/3)*p));
        end
    end
    identity = [identity;identaux(triup==1)];
    jcantor  = [jcantor;jcantaux(triup==1)];

    if g == 1
        xplotf = [0*xplot;xplot];
        plotterseq(g,1)=0;
        plotterjuk(g,1)=0;
    else
        xplotf = [xplotf;g*xplot];
    end
    
    %% Graph
    if g<=k
%         plot(xplotf,identity,'rx',xplotf,jcantor,'bo','MarkerSize',8);
%         title('Similar evolution between population individuals');
%         axis([0 ngen+1 0 2.5]);
%         xlabel('Generations');
%         ylabel('Similarity / Difference (%)');
        
        plotterseq(g,1)=mean2(identity);
        plotterjuk(g,1)=real(-(3/4)*log(1-(4/3)*plotterseq(g,1)));
    end
end;

figure
plot(plotterseq(:,1),'r');
hold on
plot(plotterjuk(:,1),'b');
xlabel('Generations');
ylabel('Distance (Seq. Averaged)');

%% Print Sequences on the terminal

fprintf('>Initial Sequence:\n');
fprintf('%s',bases{seq});
for t=1:popsize
    fprintf('\n>individuo %d\n',t);
    fprintf('%s',bases{population(t,:)});
end
fprintf('\n');

%% Save sequences on a FASTA file

% Final Population
for ff = 1:popsize
    auxseq = '';
    for ss = 1:seqsize
        auxseq = strcat(auxseq,char(bases(population(ff,ss))));
    end
    data(ff).Header   = sprintf('Sequence_%d',ff);
    data(ff).Sequence = sprintf('%s',auxseq);
end

filename = input('\n\nInsert the name of the FASTA file Exemple: ''laboratorio#5.txt''.\n');
fastawrite(filename,data);

end
