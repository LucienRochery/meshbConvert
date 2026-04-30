Build with:
```
mkdir build
cd build 
cmake .. 
make -j4
```

Run with:
```
ugrid2meshb toto.meshb toto.lb8.ugrid 
ugrid2meshb toto.lb8.ugrid toto.meshb
```
This should output the sibling `.mapbc`.
Same thing for stl:
```
stl2meshb toto.meshb toto.stl
stl2meshb toto.stl toto.meshb
```