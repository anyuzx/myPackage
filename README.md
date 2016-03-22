This python package is for my own daily scientific computation use. Some module is compiled using `Cython`. Git will ignore all the compiled `.so`, `.o` files on specific machines.

## Instruction

---

#### Compile

In order to use certain tools in this package, `Cython` code need to be compiled. Use the following command to compile the codes.

```bash
python setup.py build_ext --inplace
```

This command compile the source code and put the compiled `.so` files in the corresponding locations.

#### Known Issues

package `matplotlib` needed to be imported after this package, otherwise some cython module in this package won't link the dynamical
library correctly. For instance, if you do the following:

```python
import matplotlib.pyplot as plt
import myPackage as mp

...
```

It will show ERRORS something like this

```
...
ImportError: /usr/local/gcc/lib64/libstdc++.so.6: version 'GLIBCXX_3.4.18' not found
...
```

The reason of this bug is not known so far. Just do `import myPackage` before `import matplotlib`
