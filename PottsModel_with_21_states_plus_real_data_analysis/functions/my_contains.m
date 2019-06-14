function match_flag=my_contains(str,pattern)
L=length(str);
L1=length(pattern);
match_flag=0;

i=1;
while str(i)~='_'
    i=i+1;
end

if L-i>=L1 && sum(str(i+1:i+L1)==pattern)==L1
    match_flag=1;
end

