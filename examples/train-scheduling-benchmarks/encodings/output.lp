%%% Output
#show.
#show fahrweg(SI,F) : si(fahrweg,SI,F).
#show route(SI,ID,ID') : route(SI,(V,V')), route(SI,(V',V'')), 
                         edge(SI,ID,V,V'), edge(SI,ID',V',V'').
#show hash/1.
#show topologischeReihenfolge/3.
#show abschnittskennzeichen(SI,ID,K) : si(abschnittVorgabe,SI,K), si(fahrweg,SI,F), abschnitt(abschnittskennzeichen,F,ID,K).
#show abschnitt(folgenid,SI,ID,Name) : abschnitt(folgenid,F,ID,Name), si(fahrweg,SI,F).
#show edge/4.
