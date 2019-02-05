start((t1;t2),1).    start((t1;t2),2).    start((t1;t2),3).
end((t1;t2),11).     end((t1;t2),12).     end((t1;t2),13).
edge((t1;t2),1,4).   edge((t1;t2),2,4).   edge((t1;t2),3,4).  
edge((t1;t2),4,5).   edge((t1;t2),5,6).   edge((t1;t2),5,7).   
edge((t1;t2),6,8).   edge((t1;t2),7,9).   edge((t1;t2),8,10). 
edge((t1;t2),9,10).  edge((t1;t2),10,11). edge((t1;t2),10,12). 
edge((t1;t2),10,13).
resource(sw1,(1,4)).   resource(sw1,(2,4)).   resource(sw1,(3,4)).
resource(sw1,(4,5)).
resource(sw2,(8,10)).  resource(sw2,(9,10)).  resource(sw2,(10,11)).
resource(sw2,(10,12)). resource(sw2,(10,13)).
resource(r(V,V'),(V,V')) :- edge(_,V,V').
b(R,1)                   :- resource(R,_).
m((V,V'),1)              :- edge(_,V,V').
train(t1). train(t2).
w(T,(V,V'),0)            :- edge(T,V,V').
e(t1,(1;2;3),0).     l(t1,(1;2;3),6).     
e(t2,(1;2;3),4).     l(t2,(1;2;3),9).
e(t1,4,2).           l(t1,4,7).           
e(t2,4,5).           l(t2,4,10).
e(t1,5,3).           l(t1,5,8).           
e(t2,5,6).           l(t2,5,11).
e(t1,(6;7),4).       l(t1,(6;7),9).       
e(t2,(6;7),7).       l(t2,(6;7),12).
e(t1,(8;9),5).       l(t1,(8;9),10).      
e(t2,(8;9),8).       l(t2,(8;9),13).
e(t1,10,6).          l(t1,10,11).         
e(t2,10,9).          l(t2,10,14).
e(t1,(11;12;13),7).  l(t1,(11;12;13),13). 
e(t2,(11;12;13),10). l(t2,(11;12;13),15).
connection(t2,(6,8),t1,(7,9),1).
connection(t2,(7,9),t1,(6,8),1).
connection(t2,(6,8),t1,(6,8),1).
connection(t2,(7,9),t1,(7,9),1).

potlate(t1,(1;2;3),5,1).     potlate(t2,(1;2;3),8,1).   
potlate(t1,4,6,1).           potlate(t2,4,9,1).    
potlate(t1,5,7,1).           potlate(t2,5,10,1).
potlate(t1,(6;7),8,1).       potlate(t2,(6;7),11,1).
potlate(t1,(8;9),9,1).       potlate(t2,(8;9),12,1).
potlate(t1,10,10,1).         potlate(t2,10,13,1).
potlate(t1,(11;12;13),11,1). potlate(t2,(11;12;13),14).

penalty((2,4),1).            penalty((10,12),1).
