function plot_mode_correlation(index,NormVector)

Vector_info=size(NormVector);
N=Vector_info(1)/3;


temp_coord=zeros(N,3);
temp=NormVector(:,index)';
    
for i=1:N
temp_coord(i,1)=temp(3*i-2);
temp_coord(i,2)=temp(3*i-1);
temp_coord(i,3)=temp(3*i);
end

correlation=zeros(N,N);

for i=1:N
    for j=1:N
       correlation(i,j)=sum(temp_coord(i,:).*temp_coord(j,:)); 
    end
end

meshsize=4;
coarse_graining_corr(correlation,meshsize);
