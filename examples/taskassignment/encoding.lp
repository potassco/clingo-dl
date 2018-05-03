task(T):-executionTime(T,_,_).
machine(M):-executionTime(_,M,_).

%get maximum timepoint
maxtime(T,T1):-task(T), executionTime(T,R1,T1), 0{executionTime(T,R2,T2):R2!=R1, T2>T1}0. %usually bad encoding, fine in this case since domain predicates only
time(0..Max-1):-Max=#sum{T1,T,R1:executionTime(T,R1,T1)}. %see comment above

% decide which operation first
{ before(T1,ST1,T2,ST2) } :- assign(T1,ST1,M), assign(T2,ST2,M), (T1,ST1)!=(T2,ST2).
before(T2,ST2,T1,ST1) :- assign(T1,ST1,M), assign(T2,ST2,M), (T1,ST1)!=(T2,ST2), not before(T1,ST1,T2,ST2).  

subtasks(T,N) :- task(T), N={assign(T,ST,M)}.
1{ sequence(T,ST,O) : O=1..N,subtasks(T,N) }1 :- assign(T,ST,M).
:-not 1{sequence(T,ST,O)}1, O=1..N,subtasks(T,N).

seq(T,ST,0) :- sequence(T,ST,1).

seq(T,ST2,S+Time) :- sequence(T,ST1,N), sequence(T,ST2,N+1), executionTime(T,ST1,Time), seq(T,ST1,S), time(S+Time).
seq(T2,ST2,S+Time) :- before(T1,ST1,T2,ST2), executionTime(T1,ST1,Time), seq(T1,ST1,S), time(S+Time).
:- before(T1,ST1,T2,ST2), executionTime(T1,ST1,Time), seq(T1,ST1,S), not time(S+Time).

aux(T,ST,S)   :- seq(T,ST,S).
aux(T,ST,S-1) :- aux(T,ST,S), time(S-1).     

bind(T,ST,S) :- aux(T,ST,S), not aux(T,ST,S+1).

:- bind(T,M,S), S>B, bound(B).

#show seq/3.