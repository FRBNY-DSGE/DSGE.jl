"""
```
default_policy()
```

Creates an instance of an AltPolicy for the default policy rule
"""
function default_policy()
    AltPolicy(:default_policy, eqcond, solve, forecast_init = identity)
end
