# HEAT / Jacobi Method

This is a parallel MPI implementation of the iterative [Jacobi method](https://en.wikipedia.org/wiki/Jacobi_method) for solving linear systems of equations.

To compile this application you'll need an MPI distribution with ULFM support enabled.

To compile the SCR version of this application you'll need to install SCR (SCR 3.0.1 was used on the development of this version)

## Running the versions

After compile with `make`, Examples of how to run the different versions are below.

### NOFT

    mpirun -np NP jacobi_noft -p NR -q NC -NB QC -MB QR

### ULFM

    mpirun --with-ft=ulfm --oversubscribe -np NP jacobi_ulfm -p NR -q NC -NB QC -MB QR

### SCR

    mpirun -np NP jacobi_scr -p NR -q NC -NB QC -MB QR [--debug] [--use-scr-need-checkpoint]
    Flags:
        --debug: Print debug information (not only SCR debug default messages)
        --use-scr-need-checkpoint: Use SCR_Need_checkpoint to checkpoint the application

### Arguments

    NP: Number of processes used by the application
    NR: Number of processes per row
    NC: Number of processes per column
    QC: Number of columns
    QR: Number of rows

## Troubleshooting

If `make jacobi_scr` isn't working, change the content of the `SCRDIR` variable to the path of your SCR installation
