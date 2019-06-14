function input_data_cp=scramble_sequences(input_data)

[M,L]=size(input_data);
input_data_cp=input_data;
for k=1:L
    for j=1:floor(M/2)  % statistically,  randonly swith M/2 times have almost shift all the sequences
         l1=floor(M*rand(1,1))+1;
         l2=floor(M*rand(1,1))+1;
         input_data_cp(l1,k)=input_data(l2,k); 
         input_data_cp(l2,k)=input_data(l1,k); 
    end
end
