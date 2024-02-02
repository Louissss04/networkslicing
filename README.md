# Optimal QoS-Aware Network Slicing for Service-Oriented Networks with Flexible Routing

The repository contains the Julia code corresponding to the paper:

   Optimal QoS-Aware Network Slicing for Service-Oriented Networks with Flexible Routing

To run the code, we need to install: 
(i) GUROBI: [GUROBI Installation](https://www.gurobi.com)

## Description of the Files

1. **ReadData.jl**
   - Purpose: Used to read and store information for generating graphs.
   - Called by: Used by `test.jl`, `compare.jl`, and `CG_main.jl`.

2. **constraint.jl**
   - Purpose: Groups constraints used in optimization models across the project into functions.
   - Called by: Used by `test.jl`, `compare.jl`, and `CG_main.jl`.

3. **masterproblem.jl**
   - Purpose: Constructs the column generation master problem model.
   - Called by: Used by `CG_main.jl`.

4. **subproblem.jl**
   - Purpose: Constructs the column generation subproblem model.
   - Called by: Used by `CG_main.jl`.

5. **test.jl**
   - Purpose: Test entry point for freely testing MINLP and MILP models (switching between them, adding/removing QoS constraints).

6. **compare.jl**
   - Purpose: Test entry point for comparing the equivalence of MINLP and MILP.

7. **CG_main.jl**
   - Purpose: Test entry point for implementing the column generation algorithm.

8. **graph_fish.m**
   - Purpose: Generate test graph instances.

## Usage

1. Install dependencies:
    ```bash
    # Run the installation command
    ```

2. Run tests:
    ```bash
    # julia test.jl ntests K capacityparam
    julia test.jl 1 10 40
    ```

3. Run comparisons:
    ```bash
    julia compare.jl
    ```

4. Run the column generation algorithm:
    ```bash
    # julia CG_main.jl ntests K capacityparam
    julia CG_main.jl 1 10 40
    ```

## Notes

- Ensure Julia environment is installed and required dependencies are installed.
- Please modify the paths accordingly based on the storage location.
- Before running the column generation algorithm, it is recommended to run `test.jl` and `compare.jl` to ensure everything is functioning correctly.

## Contribution

Contributions, bug reports, and suggestions are welcome. Please email 2207740560@qq.com, Louis. 

## License

This project is licensed under the [MIT License](LICENSE).

