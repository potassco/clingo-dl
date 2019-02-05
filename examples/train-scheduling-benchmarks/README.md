# Train scheduling with clingoDL: Benchmarks
## REQUIREMENTS

- clingo-dl 1.0.0, avilable [here](https://github.com/potassco/)
- clingo 5.3.0, avilable [here](https://github.com/potassco/)
- python >= 3.6

## CONTENT

- encodings: Folder containing all ASPmDL encodings including optimization and heuristics
- instances:
	- json: Instances in json format
	- asp: Instances in ASP facts
- utils: Folder with solution converter and checker
- solve_and_check: Script running and validating one instance

## USAGE

`solve_and_check instance <config> <time limit> <max delay>`

- `instance`: ASP instance
- `config`: Configuration of clingo-dl (optional)
- `time limit`: Time limit for clingo-dl (optional)
- `max delay`: Maximum delay that is allowed for each train at a node (optional)

Example call:
`./solve_and_check.sh instances/asp/03_FWA_0.125.lp`