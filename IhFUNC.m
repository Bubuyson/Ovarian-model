function val = IhFUNC(t, dIh)
    val = 300 + 1330*exp(-(t-7- dIh)^2/19);
end