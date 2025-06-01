# checking equivalence of bezier curve using SMT Z3

Z3 is needed so you need to install it first.

Checking can be done in smt-lib format or in python script. Examples for both are provided.

directly using smt-lib format
```
z3 -smt2 bezier_curves.smt2
```

using the python z3 API
```
python3 test/snapshots/bezier_curves_z3.py
```

