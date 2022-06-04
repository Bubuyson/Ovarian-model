function val = relLH(t, RPLH)
    val = (3*(1+0.024*P4FUNC(t, 0))*RPLH)/(1+0.008*E2FUNC(t,0));
end