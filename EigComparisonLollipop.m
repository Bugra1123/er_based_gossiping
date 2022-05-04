LolComparison={"n","Lambda ER","Lambda Met","Lambda Uni","Lambda Opt"}
for n=10:10:1000
    [l_ER,l_met,l_uni,l_opt]=EigCompLol(n);
    LolComparison=[LolComparison;{n,l_ER,l_met,l_uni,l_opt}];
end