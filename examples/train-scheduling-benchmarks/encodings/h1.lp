range(T,(V,V'),E,L) :- edge(T,V,V'),e(T,V,E), l(T,V',L).
#heuristic seq(T,E,T',E') : shared(T,E,T',E'), 
                            range(T,E,L,U), 
                            range(T',E',L',U'). [L'-L - (U-U'),sign]


