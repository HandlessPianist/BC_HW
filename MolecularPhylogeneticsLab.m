%Computational Biology 2015/2016 - Molecular Phylogenetics Lab
%David Ribeiro, 72666
%Andre Pombeiro, 73268
%Claudia Vasconcelos, 73417
%% - Data insertion by the user to create the initial population
close all;clc;
disp('Laboratory #5 - Computational Biology');
option=0;

while option~=3

    fprintf('\n\nChoose one of the following options:');
    option=input('\n1. Create population.\n2. Evolve population.\n3. Exit\n');
    switch option

        case 1

            PopSize=input('\n 1.1. Enter the population size \n ');
            SeqLength=input('\n 1.2. Enter the sequence length \n');

            CreatePop(PopSize,SeqLength);

            disp('The program is paused. Press enter to continue.');
            pause;

        case 2

            PopFile=input('\n 1.1. Enter the path to file.\n ');
            NGen=input('\n 1.1. Enter the number of generations for the evolution step.\n ');
            MR=input('\n 1.1. Enter the mutation rate.\n ');
            RR=input('\n 1.1. Enter the recombination rate.\n ');
            RFL=input('\n 1.1. Enter the recombination fragment length.\n ');

            EvolvePop(PopFile,NGen,MR,RR,RFL);
            disp('The program is paused. Press enter to continue.')
            pause;

    end
end
