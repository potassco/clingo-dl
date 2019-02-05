1 { visit(T,V)      : start(T,V)    } 1 :- train(T), 2 { start(SI,V) }.
visit(T,V) :- train(T), start(T,V), { start(T,V') } 1.
1 { route(T,(V,V')) : edge(T,V,V')  } 1 :- visit(T,V), not end(T,V), 2 { edge(T,V,V') }.
route(T,(V,V')) :- visit(T,V), edge(T,V,V'), { edge(T,V,V'') } 1.
visit(T,V')     :- route(T,(V,V')).
route(T,(V,V')) :- bridge(T,ID,"True"), edge(T,ID,V,V').
      
shared(T,(V,V'),T',(U,U')) :- edge(T,V,V'), edge(T',U,U'), T!=T',
                              M = #max{ B : resource(R,(V,V')), resource(R,(U,U')), b(R,B) },
                              1 #sum{ 1,R : resource(R,(V,V')), resource(R,(U,U')) },                             
                              e(T,V,E), l(T,V',L), e(T',U,E'),
                              E <= E', E' < L+M, (E,T)<(E',T').
shared(T',(U,U'),T,(V,V')) :- shared(T,(V,V'),T',(U,U')), e(T,V,E), e(T',U,E'), (E,T)<(E',T').                        
conflict(T,E,T',E') :- shared(T,E,T',E'), T < T',
                       route(T,E), route(T',E').
{ seq(T,E,T',E') } :-  conflict(T,E,T',E').
  seq(T',E',T,E)   :-  conflict(T,E,T',E'), not seq(T,E,T',E').
:- seq(T,(V,V'),T',(U,U')),  
   route(T,(V',V'')), route(T',(U',U'')), 
   shared(T,(V',V''),T',(U',U'')), not seq(T,(V',V''),T',(U',U'')).
:- seq(T,(V,V'),T',E),   route(T,(V',V'')),
   shared(T,(V',V''),T',E), not seq(T,(V',V''),T',E).
:- seq(T,E,T',(U,U')), route(T',(U',U'')),
   shared(T,E,T',(U',U'')), not seq(T,E,T',(U',U'')).
    
&diff{ 0-(T,V) } <= -E    :- e(T,V,E), visit(T,V), relevant(e,T,V).
&diff{ (T,V)-0 } <=  L    :- l(T,V,L), visit(T,V), relevant(l,T,V).
&diff{ (T,V)-(T,V') }  <= -D :- route(T,(V,V')), E = (V,V'),
                                D=#sum{ M,m : m(E,M); W,w : w(T,E,W) }.
&diff{ (T,V')-(T',U) } <= -M :- seq(T,(V,V'),T',(U,U')), 
                                M = #max{ B : resource(R,(V,V')), resource(R,(U,U')), b(R,B) }.
&diff{ (T,V)-(T',U') } <= -W :- connection(T,(V,V'),T',(U,U'),W), route(T,(V,V')), route(T',(U,U')).
