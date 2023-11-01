A Project started in the Algorithm Engineering Course SUMMER TERM 2020 at TU-Berlin
	Luka Stärk  

# Vertex Cover vc_solver

## Overview
Vertex Cover vc_solver is a tool for solving vertex cover problems using branch and cut, reductions, and preprocessing techniques. It allows users to enable or disable different approaches, including greedy algorithms, local search, and preprocessing, to optimize the solution process.

## Installation
To use the Vertex Cover vc_solver, clone the repository or download the source code. Ensure you have Python installed on your system.

## Usage
Set the stack size limit to 64MB, install dependencies and compiles Cython files by running the script: `. ./config.sh`
(IMPORTANT: execute config inside your bash with exactly `. ./config.sh` and NOT `./config.sh`)

Run the program from the command line using the following syntax:

```bash
python main.py < input_file > output_file [options]
```

### Options
- `-h`, `--help`: Show the help message and exit.
- `-p`, `--preprocessing`: Enable or disable preprocessing (default: disabled).
- `-g`, `--greedy`: Enable or disable the greedy algorithm (default: disabled).
- `-ls`, `--local_search`: Enable or disable local search (default: disabled).
- `-lb`, `--print_lower_bound`: Enable or disable printing of lower bound (default: disabled).
- `-il`, `--ignore_limit`: Ignore limit flag (default: disabled).

## Reduction Rules
The program also allows for fine-tuning of reduction rules. Each rule can be enabled or disabled, and its frequency of application can be set. The default settings for the rules are as follows:

- **deg_1**: Enabled, Frequency: 1
- **deg_2**: Enabled, Frequency: 1
- **high_deg**: Enabled, Frequency: 1
- **buss**: Enabled, Frequency: 1
- **dom**: Enabled, Frequency: 1
- **crown**: Disabled, Frequency: 1
- **lprule**: Disabled, Frequency: 1
- **unconfined**: Disabled, Frequency: 1
- **deg3_ind_set**: Disabled, Frequency: 1

## Customization
You can customize the reduction rules and their frequencies by modifying the relevant section in the code.

