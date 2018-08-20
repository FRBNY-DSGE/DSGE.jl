Hello fello citizen of the RAN!

I have not changed the measurement and observables file spare changing
the name of the model. In bondlabor.jl, beyond changing the name, I
also altered the parameters whose values were given in the original
Julia file from Sushant and Keshav. I commented out the ones whose
values were left unspecified for easy signallying of what would need
to be changed. Then, I took \beta from being a parameter to being
a setting, and removed R_* from being a setting. I assumed we'd need
to do this because we're solving for \beta, and R_* is being held at
a fixed 1.04. Obviously, I wasn't hyper-vigilant with any of this, so
checking for typos would be wise.

Best,
Reca
8/14/2018
