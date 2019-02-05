trans(SI,ID'',(F,ID,alt)) :- jsedge(SI,ID,ID''), 2{ jsedge(SI,ID,ID') }, si(fahrweg,SI,F).
trans(SI,ID'',(F,ID''))   :- jsedge(SI,ID,ID''), { jsedge(SI,ID,ID') }1, si(fahrweg,SI,F).
trans(SI,ID,(F,ID))       :- jsedge(SI,ID,ID'), not jsedge(SI,_,ID), si(fahrweg,SI,F).
train(T) :- si(T).

edge(SI,ID,(F,ID),(F,ID'))        :- jsedge(SI,ID,ID'), { jsedge(SI,ID,ID'') }1, trans(SI,ID,(F,ID)), si(fahrweg,SI,F).
edge(SI,ID,(F,ID),(F,ID,alt))     :- jsedge(SI,ID,_), 2{ jsedge(SI,ID,ID') }, si(fahrweg,SI,F).
edge(SI,ID',(F,ID,alt),TID)       :- edge(SI,_,(F,ID,alt)), jsedge(SI,ID,ID'), jsedge(SI,ID',ID''), trans(SI,ID'',TID), si(fahrweg,SI,F). 
edge(SI,ID,(F,ID),(F,ID,end))     :- jsend(SI,ID), si(fahrweg,SI,F).
edge(SI,V,V') :- edge(SI,_,V,V').

start(T,V) :- edge(T,V,V'), not edge(T,_,V).
end(T,V')  :- edge(T,V,V'), not edge(T,V',_).

connection(SI,(V,V'),SI',(U,U'),Time) :- anschluss(aufZugfahrt,SI,K,N,SI'),
                                         anschluss(aufAbschnittskennzeichen,SI,K,N,K'),
                                         anschluss(minAnschlusszeit,SI,K,N,Time),
                                         si(fahrweg,SI,F), abschnitt(abschnittskennzeichen,F,ID,K),
                                         si(fahrweg,SI',F'), abschnitt(abschnittskennzeichen,F',ID',K'),
                                         edge(SI,ID,V,V'), edge(SI',ID',U,U').   

delta(B)   :- maxBandabweichung(B).
delta(mbw) :- not maxBandabweichung(_). 
ausMax'(SI,ID,Max+MBW) :- ausMax(SI,ID,Max), delta(MBW).
einMax'(SI,ID,Max+MBW) :- einMax(SI,ID,Max), delta(MBW).

relevant(e,SI,V)  :-  jsstart(SI,ID), trans(SI,ID,V).
relevant(e,SI,V)  :-  jsend(SI,ID), trans(SI,ID,V).
relevant(e,SI,V)  :-  abschnittVorgabe(einMin,SI,K,_), 
                          si(fahrweg,SI,F), 
                          abschnitt(abschnittskennzeichen,F,ID,K),
                          trans(SI,ID,V).
relevant(l,SI,V') :-  abschnittVorgabe(ausMax,SI,K,_),
                      si(fahrweg,SI,F),
                      abschnitt(abschnittskennzeichen,F,ID,K),
                      edge(SI,ID,V,V').
relevant(l,SI,V)  :-  si(fahrweg,SI,F), trans(SI,ID,V),
                      abschnittVorgabe(einMax,SI,K,_),
                      abschnitt(abschnittskennzeichen,F,ID,K).
relevant(e,SI,V') :-  edge(SI,ID,V,V'), si(fahrweg,SI,F),
                      abschnitt(abschnittskennzeichen,F,ID,K),
                      abschnittVorgabe(ausMin,SI,K,_).

e(SI,V,E)      :- einMin(SI,ID,E), trans(SI,ID,V).
l(SI,V,L)      :- einMax'(SI,ID,L), trans(SI,ID,V).
l(SI,V',L)     :- ausMax'(SI,ID,L), edge(SI,ID,V,V').
e(SI,V',E)     :- ausMin(SI,ID,E), edge(SI,ID,V,V').

siid(SI,(ID;ID'))     :-    jsedge(SI,ID,ID').
m((V,V'),FZ) :- siid(SI,ID), si(fahrweg,SI,F), abschnitt(minimaleFahrzeit,F,ID,FZ),
                edge(SI,ID,V,V').
w(SI,(V,V'),HZ) :- siid(SI,ID), si(fahrweg,SI,F), abschnitt(minimaleFahrzeit,F,ID,FZ),
                   edge(SI,ID,V,V'), 
                   abschnittVorgabe(minHaltezeit,SI,K,HZ),
                   abschnitt(abschnittskennzeichen,F,ID,K).
                               
resource(R,(V,V')) :- ressourcenbelegung(ressource,F,ID,R),
                      edge(SI,ID,V,V'),
                      si(fahrweg,SI,F).
b(R,B) :- ressource(freigabezeit,R,B).

penalty((V,V'),P) :- abschnitt(penalty,F,ID,P), si(fahrweg,SI,F),
                     edge(SI,ID,V,V').
                     
potlate(SI,V,M+T,W*U) :- abschnittVorgabe(einVerspaetungsfaktor,SI,K,W), si(fahrweg,SI,F),
                         abschnitt(abschnittskennzeichen,F,ID,K), abschnittVorgabe(einMax,SI,K,M), dpenalty(T,U),
                         edge(SI,ID,V,V').

potlate(SI,V',M+T,W*U) :- abschnittVorgabe(ausVerspaetungsfaktor,SI,K,W), si(fahrweg,SI,F),
                          abschnitt(abschnittskennzeichen,F,ID,K), abschnittVorgabe(ausMax,SI,K,M), dpenalty(T,U),
                          edge(SI,ID,V,V').
