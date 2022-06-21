# Learn libIgl
## build igl environment on win10
- git clone libigl, put it at `~/work_path/libigl`

    ```bash
    git clone https://github.com/libigl/libigl.git
    ```

- set this repo path to `~/work_path/example-project`

- build

    ```bash
    mkdir build
    cd build
    cmake "Visual Studio 15 2017" -A x64 ..
    cmake --build . --config Release -- /maxcpucount:8
    ```