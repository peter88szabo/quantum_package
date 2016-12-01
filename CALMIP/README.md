Quantum Package CALMIP benchmark
================================

Introduction
------------

The Quantum Package is a program for ab-initio quanum chemistry which is
designed to use both supercomputers and grid/cloud environments.  It is
parallelized with OpenMP for single-node parallelism, and ZeroMQ for both
single-node and multi-node parallelism.

Parallelism is expressed as tasks that can be started in new threads to benefit
from shared memory, or by new independent processes that can have multiple
threads. Here, we don't use the SPMD paradigm but the MPMD. Typically, a large
parallel run is started as one single multi-threaded porcess per node which can
be started by multiple SLURM sumbissions. The first job starts the task server
and spawns the master process, and all the other jobs start slave processes
that connect to the task server and execute tasks, sending the results to the
master process. Slaves can be added/removed dynamically during the run.

The computational processes are written in Fortran using the IRPF90 code generator.
The task server is written in OCaml. All the processes communicate with ZeroMQ
(asynchronous TCP).


Compiling the code
------------------

Download the code here : https://codeload.github.com/scemama/quantum_package/zip/calmip`
or using git:

```bash
git clone -b calmip https://github.com/scemama/quantum_package.git
```

and go into the ``quantum_package`` directory.

To compile the code, you will need a configuration file that fits your
architecture.  Examples are in the `config` directory. Copy a file which fits
your architecture, modify it and use it to conifgure:

```bash
cd config
cp ifort.cfg calmip.cfg
vim calmip.cfg             # Edit the config file
cd ..
./configure --production config/calmip.cfg
```

During the configuration, the missing dependencies are downloaded and compiled.
Before compiling or running the program, it is necessary to have some environment
variables set. To do this :

```bash
source quantum_package.rc
```

Now, install the modules necessary for the benchmark:

```bash
cd $QP_ROOT
qp_module.py install Hartree_Fock Full_CI_ZMQ MP2
```

Compile the program

```bash
cd $QP_ROOT
ninja
```

Benchmark
---------

Go into the `CALMIP` directory and extract the data set :

```bash
cd ${QP_ROOT}/CALMIP
tar -zxvf FeO4.tar.gz
```

Run the benchmark :

```
qp_run fci_zmq FeO4 > FeO4.out
```

There is a reference file `FeO4.ref` in the directory to check that the results are correct.
The maximum required memory s 96 GB and the run requires roughly 40 minutes on 64 cores.


If there is a problem with the dataset, it can be re-generated using the script `generate_dataset.sh`.
