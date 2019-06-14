function  WW_seq=fasta_to_sequence

%% mapping table
number2letter_map(1) = '-';
number2letter_map(2) = 'A';
number2letter_map(3) = 'C';
number2letter_map(4) = 'D';
number2letter_map(5) = 'E';
number2letter_map(6) = 'F';
number2letter_map(7) = 'G';
number2letter_map(8) = 'H';
number2letter_map(9) = 'I';
number2letter_map(10) = 'K';
number2letter_map(11) = 'L';
number2letter_map(12) = 'M';
number2letter_map(13) = 'N';
number2letter_map(14) = 'P';
number2letter_map(15) = 'Q';
number2letter_map(16) = 'R';
number2letter_map(17) = 'S';
number2letter_map(18) = 'T';
number2letter_map(19) = 'V';
number2letter_map(20) = 'W';
number2letter_map(21) = 'Y';

% 
% M=240; L=92; amino_N=21;
% WW_seq=zeros(M,L);
% for j=1:M
%     temp=algn(j,:);
%     for k=1:L
%         l=1;
%         while temp(k)~=number2letter_map(l) && l<amino_N
%             l=l+1;
%         end
%         
%         if temp(k)==number2letter_map(l)
%             WW_seq(j,k)=l;
%         else 
%             WW_seq(j,k)=1;
%         end
%         
%     end
%     
% end


%%
[WW_head,WW_seq_alphabeta]=fastaread('PF00397_uniprot.fasta');

%%
amino_N=21;
M=length(WW_seq_alphabeta);
L=length(WW_seq_alphabeta{1});
WW_seq=zeros(M,L);

for j=1:M
    temp=WW_seq_alphabeta{j};
    for k=1:L
        l=1;
        while temp(k)~=number2letter_map(l) && l<amino_N
            l=l+1;
        end
        
        if temp(k)==number2letter_map(l)
            WW_seq(j,k)=l;
        else 
            WW_seq(j,k)=1;
        end
        
    end
    
end



