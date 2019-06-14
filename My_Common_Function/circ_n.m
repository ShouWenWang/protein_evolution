    function m=circ_n(n,N)
        m=mod(n,N);
        m(mod(n,N)==0)=N;
    end