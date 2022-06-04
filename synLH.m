function val = synLH(t, dE, dP)
    val = (1400 + (95900*E2FUNC(t, dE)^8/(360^8+E2FUNC(t, dE)^8)))/(1+ P4FUNC(t, dP)/26);
end