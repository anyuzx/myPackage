This python package is for my own daily scientific computation use. Some module is compiled using `Cython`. Git will ignore all the compiled `.so`, `.o` files on specific machines.

## Instruction

---

#### Compile

In order to use certain tools in this package, `Cython` code need to be compiled. Use the following command to compile the codes.

```bash
python setup.py build_ext --inplace
```

This command compile the source code and put the compiled `.so` files in the corresponding locations.

