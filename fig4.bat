@java -Xmx1024m -jar doubleModel.jar -j 1 33 -rp 100000 -pdflim 1e-3 1e3 --limits 1e-3 0.999 -op 300 -k 0.03 -ot 0. -pars 0.1 3 3 300 -m 0
@java -Xmx1024m -jar doubleModel.jar -j 1 33 -rp 100000 -pdflim 1e-3 1e3 --limits 1e-3 0.999 -op 300 -k 0.03 -ot 1. -pars 0.1 3 3 300 -m 1
@java -Xmx1024m -jar doubleModel.jar -j 1 99 -rp 100000 -pdflim 1e-3 1e3 --limits 1e-3 0.999 -op 300 -k 0.03 -ot 2. -pars 0.1 3 3 300 -m 2
@java -Xmx1024m -jar doubleModel.jar -j 1 333 -rp 100000 -pdflim 1e-3 1e3 --limits 1e-3 0.999 -op 300 -k 0.03 -ot 3. -pars 0.1 3 3 300 -m 4
@java -Xmx1024m -jar doubleModel.jar -j 1 999 -rp 100000 -pdflim 1e-3 1e3 --limits 1e-3 0.999 -op 300 -k 0.03 -ot 4. -pars 0.1 3 3 300 -m 8