function val = E2FUNC(t, dE)
    val = 300 - 240*(t+1-dE)^2/(3+(t+1-dE)^2)+90*exp(-(t-8-dE)^2/10);
end