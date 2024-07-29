* Added
- CMakeLists

Dependencies:
- spdlog to log
- fmt (come with spdlog)

Executables
- work_factor_computation_parallel
To spwan threads that compute the work_factor

* Structure
- include/utils
All the hpp headers
- src/utils
Al the cpp corresponding to the prev headers
- src/tools
All the output tools, that is, executables to use

* Compile

```sh
mkdir build && cd build 
cmake ..
make -j
```

** To create the binaries (inside the local `bin` directory)
Inside `build`

```sh
cmake -DCMAKE_INSTALL_PREFIX=.. ..
make -j
make install
```

Then, you can execute the files as ./bin/<tool_name>

* TODOs
The tools should not hardcode their paths (see f.e. work_factor_computation_parallel)
