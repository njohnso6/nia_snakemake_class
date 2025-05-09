# Alternative Actions

## Run

Shell isn't the only directive available for actions. Pure Python may be used if desired.

In order to do so, the `run` directive is needed.

```snakemake
rule example_run:
    input:
        "data/input.txt"
    output:
        "data/output.txt"
    run:
        # You can write Python code here
        with open(input[0], 'r') as infile:
            data = infile.read()
        
        # Modify the data somehow
        modified_data = data.upper()

        # Save the result to the output file
        with open(output[0], 'w') as outfile:
            outfile.write(modified_data)
```

Notice variables such as input and output are referred to as dictionaries of content.

Generally, this is best for small tasks as a number of functions don't work when using run, such as conda and singularity environments, and logging can be messed up.

## Script

Another option is called `script`. This helps keep your workflow cleaner, especially when the logic is more complex.

```snakemake
rule example_script:
    input:
        "data/input.txt"
    output:
        "data/output.txt"
    script:
        "scripts/process_data.py"
```

In a separate script called `process_data.py` we can then put

```python
# This script receives `input`, `output`, `params`, etc. from Snakemake

with open(snakemake.input[0], 'r') as infile:
    data = infile.read()

# Process the data (e.g., converting to uppercase)
modified_data = data.upper()

# Write to output
with open(snakemake.output[0], 'w') as outfile:
    outfile.write(modified_data)
```

Notice we refer to the different parts of the rule using a snakemake object.
