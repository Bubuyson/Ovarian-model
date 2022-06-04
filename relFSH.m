function val = relFSH(t, RPFSH)

    val = 65*(1+4*P4FUNC(t,0)*RPFSH)/(1+0.007*E2FUNC(t, 0)^2);

end