function pop = CreatePop(PopSize,SeqLength)

InitPop=zeros(PopSize,SeqLength); % Initialization of a population of size PopSize
bases={'A','C','G','T'}; % Nucleotides
corr={1,2,3,4}; % Nucleotides dictionary
sMap=containers.Map(corr,bases); % Nucleotides mapping

Pop=randi([1 4],1,SeqLength); % Generation of a random sequence for the initial population
    for k=1:PopSize
        InitPop(k,:)=Pop(1,:);
    end

fasta=struct('Header',{''},'Sequence',{''});

for k=1:PopSize
fasta(k).Header=strcat('>individuo_',num2str(k));
aux='';
for l=1:SeqLength
aux=strcat(aux,sMap(InitPop(k,l)));
end
fasta(k).Sequence=aux;
end

%Saving FASTA file

filename = input('\n\nInsert the name of the FASTA file. Example: ''initial_pop_trial_x.txt''.\n');
fastawrite(filename,fasta);
end
