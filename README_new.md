* Added
- CMakeLists

Dependencies:
- spdlog to log
- fmt (come with spdlog)

Executables
- work_factor_computation
Compute the work factor. It accepts either a json or a plain set of parameters

* Structure
- include/utils
All the hpp headers
- src/utils
Al the cpp corresponding to the prev headers
- src/tools
All the output tools, that is, executables to use

* Compile
To create the binaries (inside the local `bin` directory)
Inside `build`

```sh
cmake -DCMAKE_INSTALL_PREFIX=.. ..
make install -j
```

Then, you can execute the files as ./bin/<tool_name>
