node(T,(V;V'))   :- edge(T,V,V').
node_pos(T,V,0)  :- start(T,V).
node_pos(T,V',M+1) :- node(T,V'), not start(T,V'),
                      M = #max{ P : node_pos(T,V,P), edge(T,V,V')},
                      M!=#inf.
last_node(T,M)   :- train(T), M = #max{ P : node_pos(T,V,P), end(T,V) },
                    M!=#inf.
#heuristic late(T,V,U,W) : potlate(T,V,U,W),
                           node_pos(T,V,P),
                           last_node(T,Max). [Max-P,false]
