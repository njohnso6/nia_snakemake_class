# Setting up Resources

## In rules

## Benchmarking

- dynamic resource assignment -> rule:
    input:    ...
    output:   ...
  - benchmarking
    resources:
        mem_mb=lambda wc, input: max(2.5 * input.size_mb, 300)
    shell:
        "..."

