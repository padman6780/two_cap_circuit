function Out = FunCm_drive(Cm0,USfreq, t, ZONE)
    if ZONE == 2
        Cm = -0.008 * sin(2*pi*USfreq*t);
        Cm = Cm0 + Cm;
        Out = Cm;
    else
        Out = Cm0;
    end
end